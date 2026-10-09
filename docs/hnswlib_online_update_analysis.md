# hnswlib의 HNSW 검색·온라인 갱신 알고리즘 분석

- 분석 대상: [nmslib/hnswlib](https://github.com/nmslib/hnswlib), 커밋 `ca42672` (2026-09-14)
- 핵심 파일: `hnswlib/hnswalg.h` (1,616줄), `hnswlib/visited_list_pool.h`
- 본문에서 `L123`처럼 적은 줄 번호는 위 커밋의 `hnswalg.h` 기준입니다.
- 재현 방법: `git clone https://github.com/nmslib/hnswlib third_party/hnswlib && pip install ./third_party/hnswlib` 후 `experiments/hnsw_online_update/`의 스크립트를 실행합니다.

## 0. 요약

hnswlib는 HNSW 원 논문(Malkov & Yashunin, 2018)의 저자가 관리하는 레퍼런스 구현이며, 헤더 파일만으로 구성되어 있습니다. 인덱스를 다시 만들지 않고도 다음 네 가지 온라인 연산을 지원합니다.

| 연산 | API | 그래프를 실제로 바꾸는가 |
|---|---|---|
| 삽입 | `addPoint(vec, label)` | 예. 새 노드를 만들고 양방향 간선을 연결합니다 |
| 삭제 | `markDelete(label)` | 아니요. 플래그 1비트만 켭니다 (tombstone) |
| 값 갱신 | 같은 `label`로 `addPoint` | 예. 이웃의 간선을 국소적으로 다시 연결합니다 (`updatePoint`) |
| 빈자리 재사용 | `addPoint(..., replace_deleted=true)` | 예. 삭제된 슬롯에 새 벡터를 넣고 `updatePoint`를 실행합니다 |

분석에서 확인한 핵심 사항은 다음과 같습니다.

1. **삭제는 논리적 삭제뿐입니다.** 삭제된 노드도 검색 중에 계속 탐색 경로로 쓰이고, 결과 집합에만 들어가지 않습니다. 그래서 연결성은 유지되지만, 삭제 비율이 높아지면 검색 비용이 늘어납니다. 실험에서는 삭제 비율이 90%일 때 쿼리 지연이 약 3.2배로 늘었습니다.
2. **물리적 삭제와 압축(compaction)이 없습니다.** `cur_element_count`는 줄어들지 않습니다. `replace_deleted`를 켜지 않으면 삭제된 슬롯이 차지한 메모리를 영구히 회수하지 못합니다.
3. **갱신은 국소적인 재연결입니다.** 갱신 대상 노드의 1-hop 이웃만 2-hop 후보 집합에서 이웃 목록을 다시 고르고, 대상 노드 자체는 처음 삽입할 때와 같은 방식으로 다시 연결합니다. 대상 노드를 가리키기만 하고 대상 노드는 가리키지 않는 노드(단방향 inbound 이웃)의 간선은 수정하지 않으므로 오래된 간선이 남습니다.
4. **갱신 비용은 삽입보다 큽니다.** 실험에서 `replace_deleted` 경로로 넣은 점은 신규 삽입보다 점당 약 5.5배 느렸습니다.
5. **품질 저하는 작습니다.** 전체 데이터의 110%에 해당하는 양을 교체한 뒤에도 `ef=50` 이상에서는 새로 만든 인덱스와 recall@10 차이가 0.002 이내였습니다. `ef=20`처럼 낮은 ef에서만 약 1.5%p 차이가 났습니다.
6. **동시성 보장이 제한적입니다.** 삽입끼리는 동시에 실행해도 안전하지만, 검색과 삽입·갱신을 동시에 실행하는 것은 안전하지 않습니다(README 명시). 검색 경로는 인접 리스트를 잠금 없이 읽습니다.

## 1. 자료구조와 메모리 레이아웃

### 1.1 주요 파라미터 (`L95-L154`)

| 변수 | 값 | 의미 |
|---|---|---|
| `M_` | 사용자 지정 (기본 16) | 새 노드가 만드는 간선 수, 상위 계층의 최대 차수 |
| `maxM_` | `M_` | 1계층 이상의 최대 차수 |
| `maxM0_` | `2 * M_` | 0계층의 최대 차수 |
| `ef_construction_` | `max(efc, M_)` | 삽입할 때 쓰는 후보 리스트 크기 |
| `ef_` | 기본 10 | 검색할 때 쓰는 후보 리스트 크기. 실제로는 `max(ef_, k)`를 씁니다 |
| `mult_` | `1 / ln(M_)` | 레벨을 샘플링할 때 쓰는 정규화 상수 |

### 1.2 0계층 레코드

0계층은 `max_elements * size_data_per_element_` 크기의 연속된 블록 하나에 들어 있습니다(`L132`). 노드 한 개의 레코드 구조는 다음과 같습니다.

```
offset 0                                     offsetData_          label_offset_
| cnt(2B) | flags(1B) | pad(1B) | nbr[0..2M-1] (4B each) | vector (dim*4B) | label (8B) |
            └ bit0 = DELETE_MARK
```

- 이웃 수(`cnt`)는 하위 2바이트에, 삭제 플래그는 세 번째 바이트의 최하위 비트에 들어 있습니다(`L1017-L1030`). `setListCount`는 2바이트만 덮어쓰므로 삭제 플래그를 건드리지 않습니다.
- 간선, 벡터, 라벨을 한 레코드에 붙여 두었기 때문에, 노드를 확장할 때 캐시 지역성이 좋습니다. 검색 루프는 다음 이웃의 visited 슬롯과 벡터를 `_mm_prefetch`로 미리 읽어 둡니다(`L396-L414`).
- 예: `M=16`, `dim=64` float일 때 레코드 하나는 4 + 128 + 256 + 8 = **396 B**입니다.

### 1.3 상위 계층

1계층 이상은 노드마다 `malloc(size_links_per_element_ * level + 1)`로 따로 할당합니다(`L1344`). 레벨 하나당 `4 + 4*M` 바이트를 씁니다. 레벨이 0인 노드는 이 할당을 하지 않습니다.

### 1.4 레벨 샘플링 (`L227-L231`)

```
level = floor(-ln(U(0,1)) * mult_),  mult_ = 1/ln(M)
```

이 식에서 P(level ≥ l) = M^(-l)입니다. `M=16`이면 노드의 1/16이 1계층 이상에, 1/256이 2계층 이상에 들어갑니다. 계층이 올라갈 때마다 노드 수가 1/M로 줄어드는 skip list와 같은 구조입니다.

### 1.5 방문 표시 (`visited_list_pool.h`)

- 검색 한 번마다 길이 `max_elements`인 `uint16` 배열을 풀에서 빌려 씁니다.
- 배열을 매번 0으로 지우지 않고 태그 값 `curV`를 1씩 올립니다. `visited[id] == curV`이면 이번 검색에서 이미 방문한 노드입니다. 태그가 65535를 넘어 0으로 돌아올 때만 `memset`을 실행합니다.
- 풀은 동시에 실행되는 검색 수만큼 커집니다. 스레드 하나당 `2 * max_elements` 바이트를 추가로 씁니다.

## 2. 검색 알고리즘

`searchKnnNoExceptions` (`L1415-L1472`)는 두 단계로 구성됩니다.

### 2.1 단계 1: 상위 계층 탐욕 하강 (`L1426-L1451`)

```
cur = enterpoint_node_
for level = maxlevel_ down to 1:
    repeat:
        for each nbr in links(cur, level):
            if dist(q, nbr) < dist(q, cur): cur = nbr   # 개선되면 즉시 이동
    until 이번 패스에서 개선 없음
```

- 각 계층에서 후보 리스트 크기가 1인 탐욕 검색을 합니다. 개선이 멈춘 지점이 다음 계층의 진입점입니다.
- **삭제 플래그와 필터를 확인하지 않습니다.** 삭제된 노드도 경유지로 그대로 사용합니다.
- 진입점 `enterpoint_node_`는 삭제되어도 바뀌지 않습니다.

### 2.2 단계 2: 0계층 빔 검색 `searchBaseLayerST` (`L335-L472`)

두 개의 우선순위 큐를 씁니다.

- `candidate_set`: 확장할 후보. 거리의 부호를 뒤집어 넣으므로 가장 가까운 노드가 맨 위에 옵니다.
- `top_candidates`: 결과 후보, 크기 ≤ ef. 최대 힙이므로 맨 위가 현재 결과 중 가장 먼 노드이고, 그 거리가 `lowerBound`입니다.

```
push ep into candidate_set
if ep가 유효(삭제되지 않았고 필터 통과): push ep into top; lowerBound = d(ep)
else: lowerBound = +inf
while candidate_set not empty:
    c = 가장 가까운 후보
    if d(c) > lowerBound and |top| == ef: break        # 종료 조건
    pop c
    for nbr in links0(c):
        if visited: continue
        d = dist(q, nbr)
        if |top| < ef or d < lowerBound:
            push nbr into candidate_set                 # 삭제된 노드도 확장 대상
            if nbr가 유효: push nbr into top           # 결과에는 유효한 노드만
            while |top| > ef: pop top
            lowerBound = top.max
```

#### 템플릿 분기 `bare_bone_search`

- 삭제된 노드가 하나도 없고(`num_deleted_ == 0`) 필터도 없으면 `bare_bone_search=true` 버전을 씁니다(`L1454`). 이 버전은 삭제·필터 확인을 컴파일 시점에 제거합니다.
- 삭제가 **한 건이라도** 있으면 모든 쿼리가 확인 로직이 들어간 느린 경로로 바뀝니다. 느린 경로에서는 후보마다 `isMarkedDeleted`를 확인하고, 필터가 있으면 `getExternalLabel`을 읽은 뒤 사용자 함수를 호출합니다.
- 두 버전은 종료 조건도 다릅니다. bare-bone은 `d(c) > lowerBound`만 확인하고, 느린 경로는 `|top| == ef`까지 함께 확인합니다. 삭제된 노드가 많으면 `top`이 ef개로 늦게 차기 때문에, 느린 경로는 결과를 다 채울 때까지 더 오래 탐색합니다.

#### 삭제된 노드를 통과시키는 이유와 비용

삭제된 노드를 그래프에서 빼면 그 노드를 경유하던 경로가 끊어져 도달할 수 없는 영역이 생길 수 있습니다. hnswlib는 삭제된 노드를 **경유지로는 계속 쓰고 결과에서만 제외하는** 방식으로 이 문제를 피합니다. 그 대가로 삭제 비율이 높을수록 유효한 결과 ef개를 모으기 위해 더 많은 노드를 확장해야 합니다(4.2절 실험 참고).

필터(`isIdAllowed`)도 같은 방식으로 처리합니다. 필터를 통과하지 못한 노드도 확장은 하므로, 필터의 선택도가 낮으면 탐색량이 크게 늘어납니다.

#### 사용자 정의 종료 조건 (`searchStopConditionClosest`, `L1488-L1539`)

`BaseSearchStopCondition`을 넘기면 종료 판단(`should_stop_search`), 후보 채택(`should_consider_candidate`), 결과 축소(`should_remove_extra`)를 사용자 정의 로직으로 바꿀 수 있습니다. `stop_condition.h`에는 다음 두 가지 구현이 들어 있습니다.

- `EpsilonSearchStopCondition`: 반경 ε 안에 있는 점을 찾는 range search
- `MultiVectorSearchStopCondition`: 한 문서가 여러 벡터를 가질 때 문서 단위로 상위 k개를 반환하는 검색. RAG에서 청크 단위로 색인하고 문서 단위로 검색할 때 바로 쓸 수 있습니다.

## 3. 갱신(update) 알고리즘

### 3.1 삽입 `addPointWithLevel` (`L1279-L1411`)

```
1. label_lookup_에 label이 이미 있으면 → updatePoint로 분기 (3.4절)
2. cur_c = cur_element_count++         # 내부 ID는 증가만 하는 일련번호
3. lock(link_list_locks_[cur_c])        # 삽입이 끝날 때까지 유지
4. curlevel = getRandomLevel()
5. if curlevel > maxlevel_: global 잠금을 삽입이 끝날 때까지 유지
6. 레코드 초기화, 벡터와 라벨 복사, 상위 계층 링크 메모리 할당
7. maxlevel_ … curlevel+1 계층: 2.1절과 같은 탐욕 하강(노드별 잠금 사용)
8. min(curlevel, maxlevel_) … 0 계층마다:
       top = searchBaseLayer(currObj, q, level)    # ef_construction 크기의 빔 검색
       if 진입점이 삭제됨: top에 진입점을 강제로 추가
       currObj = mutuallyConnectNewElement(...)     # 3.2절
9. curlevel > maxlevel_이면 진입점을 cur_c로 바꿉니다
```

삽입용 `searchBaseLayer` (`L245-L331`)는 검색용 `searchBaseLayerST`와 구조가 같지만 두 가지가 다릅니다.

- 확장하는 노드마다 `link_list_locks_[node]`를 잡습니다. 다른 삽입 스레드가 같은 리스트를 동시에 고칠 수 있기 때문입니다.
- 후보 리스트 크기로 `ef_construction_`을 쓰고, 0계층이 아닌 계층에서도 동작합니다.

**삭제된 진입점 처리** (`L1379-L1390`): 진입점이 삭제된 상태이면 `searchBaseLayer`의 결과에 진입점이 들어가지 않습니다. 코드는 이 경우 진입점을 후보에 강제로 넣어 새 노드가 진입점과 연결되도록 합니다. 진입점은 모든 검색이 출발하는 노드이므로 연결성을 유지하려는 목적으로 보입니다. 그 결과 새 노드의 간선 한 자리를 삭제된 노드가 차지할 수 있습니다.

### 3.2 양방향 연결 `mutuallyConnectNewElement` (`L538-L662`)

1. **정방향 간선**: `getNeighborsByHeuristic2(top, M_)`로 이웃을 최대 `M_`개 고르고 새 노드의 리스트에 씁니다. 0계층의 상한은 `2M`이지만, 새 노드는 처음에 **M개까지만** 연결합니다. 나머지 자리는 이후에 삽입되는 노드가 역방향 간선으로 채웁니다.
2. **역방향 간선**: 선택된 이웃 `n`마다 `lock(n)`을 잡고 다음을 수행합니다.
   - `n`의 리스트에 여유가 있으면 새 노드를 뒤에 덧붙입니다.
   - 리스트가 꽉 차 있으면 `{n의 기존 이웃} ∪ {새 노드}`를 후보로 다시 휴리스틱을 실행하고, `n`의 리스트를 결과로 덮어씁니다.
3. 반환값은 선택된 이웃 중 가장 가까운 노드이며, 다음 하위 계층의 진입점으로 씁니다.

2단계에는 품질상 주목할 점이 두 가지 있습니다.

- 휴리스틱이 새 노드를 탈락시키면 간선이 단방향으로 남습니다. 새 노드는 `n`을 가리키지만 `n`은 새 노드를 가리키지 않습니다.
- 휴리스틱은 상한보다 **적은** 수를 돌려줄 수 있습니다. hnswlib는 원 논문의 `keepPrunedConnections` 옵션을 구현하지 않았기 때문에, 탈락한 기존 간선은 그대로 사라집니다. 그 결과 어떤 노드가 inbound 간선을 모두 잃을 수 있습니다. `checkIntegrity()`(`L1542`)는 모든 노드의 inbound 간선 수가 0보다 큰지 assert로 확인하지만, 삽입 과정에서 이 조건을 강제하지는 않습니다.

### 3.3 이웃 선택 휴리스틱 `getNeighborsByHeuristic2` (`L475-L515`)

원 논문의 Algorithm 4(`SELECT-NEIGHBORS-HEURISTIC`)에서 `extendCandidates`와 `keepPrunedConnections`를 끈 형태입니다. HNSW 논문 이후 RNG(relative neighborhood graph) pruning이라고도 부르는 규칙입니다.

```
if |C| < M: return C                     # 후보가 M개보다 적으면 가지치기하지 않음
R = []
for c in C (기준점 q에 가까운 순서):
    if |R| >= M: break
    if 모든 r ∈ R에 대해 dist(r, c) >= dist(q, c): R.append(c)
return R
```

이미 고른 이웃 `r`이 `q`보다 `c`에 더 가까우면 `c`를 버립니다. `q → r → c` 경로로 `c`에 도달할 수 있으므로 직접 간선이 중복이라고 보는 것입니다. 이 규칙 덕분에 간선이 한 방향으로 몰리지 않고 여러 방향으로 퍼지며, 군집 사이를 잇는 장거리 간선이 남습니다. 비용은 `|R| ≤ M`에 대해 후보 하나당 최대 M번의 거리 계산이므로 O(|C|·M)입니다.

### 3.4 값 갱신 `updatePoint` (`L1109-L1185`)

이미 존재하는 `label`로 `addPoint`를 호출하면 이 함수가 실행됩니다. 라이브러리 내부 호출은 항상 `updateNeighborProbability = 1.0`을 넘깁니다(`L1085`, `L1305`).

```
1. memcpy(벡터 영역, 새 벡터)                       # 잠금 없이 덮어씀
2. for layer in 0 .. level(x):
       N1 = links(x, layer)                          # x의 1-hop 이웃
       sCand = {x} ∪ N1
       sNeigh = N1 중 확률 p로 고른 노드             # p = 1.0이면 전부
       for n in sNeigh: sCand ∪= links(n, layer)     # 2-hop까지 후보에 추가
       for n in sNeigh:
           C = sCand \ {n} 중 n에 가까운 min(efc, |sCand|)개
           links(n, layer) = Heuristic(C, layer==0 ? 2M : M)   # n의 리스트를 통째로 교체
3. repairConnectionsForUpdate(x)                     # 3.5절
```

2단계는 x가 이동했을 때 **x의 기존 이웃**들이 x를 계속 이웃으로 둘지, x 대신 다른 노드를 고를지 다시 판단하는 과정입니다. 이웃 n의 새 리스트는 n 자신의 기존 이웃과 x의 2-hop 노드 중에서 고르므로, n이 기존에 가지던 좋은 간선도 후보에 남습니다.

**계산량**: `M=16`이면 0계층에서 |N1| ≤ 32이고 |sCand| ≤ 1 + 32 + 32·32 ≈ 1,057입니다. 이웃 32개마다 후보 약 1,057개와 거리를 계산하므로 2단계만으로 거리 계산이 약 3.4만 번입니다. 여기에 휴리스틱 비용과 3단계의 `ef_construction` 빔 검색 비용이 더해집니다. 신규 삽입은 빔 검색과 역방향 간선 처리만 하므로, 갱신이 몇 배 더 비쌉니다. 4.1절 실험에서 측정한 차이는 약 5.5배였습니다.

**한계**:

- **단방향 inbound 간선이 남습니다.** y → x 간선은 있지만 x → y 간선이 없는 노드 y는 N1에 들어가지 않으므로 리스트를 다시 고르지 않습니다. y는 x가 멀리 이동한 뒤에도 x를 가리킵니다. 검색 결과가 틀리지는 않지만, 그 간선은 사실상 임의의 장거리 간선이 되어 y의 이웃 자리 하나를 차지합니다.
- **x의 레벨을 다시 뽑지 않습니다.** 갱신된 노드는 이전 레벨을 그대로 유지합니다. 레벨은 벡터 값과 무관하게 무작위로 뽑으므로 분포상 문제는 없습니다.
- **원자성이 없습니다.** 1단계에서 벡터를 잠금 없이 바꾸므로, 동시에 실행 중인 검색은 일부만 바뀐 벡터를 읽을 수 있습니다. README가 검색과 갱신을 동시에 실행하지 말라고 하는 이유 중 하나입니다.

### 3.5 재연결 `repairConnectionsForUpdate` (`L1188-L1262`)

x를 새 위치에 다시 삽입하는 과정입니다. 신규 삽입(3.1절의 7~8단계)과 거의 같습니다.

1. 진입점에서 x의 레벨보다 높은 계층까지 탐욕 하강을 합니다.
2. x의 레벨부터 0계층까지 `searchBaseLayer`를 실행하고, 결과에서 x 자신을 제외합니다(self-loop 방지).
3. `mutuallyConnectNewElement(..., isUpdate=true)`를 호출합니다. x의 리스트를 덮어쓰고, 이미 x를 가리키는 이웃에는 역방향 간선을 중복으로 추가하지 않습니다(`L606-L617`).

### 3.6 삭제 `markDelete` (`L969-L982`, `L1579-L1593`)

```
lock(label 단위 잠금)
id = label_lookup_[label]
flags(id) |= DELETE_MARK;  num_deleted_++
if allow_replace_deleted_: deleted_elements.insert(id)
```

- 간선, 벡터, 라벨 매핑을 모두 그대로 둡니다. `unmarkDelete`로 되돌릴 수 있습니다.
- 삭제된 노드는 이웃의 차수 한 자리를 계속 차지합니다. 다른 노드에 역방향 간선을 추가할 때 휴리스틱이 삭제된 노드를 탈락시키는 경우에만 자연스럽게 정리됩니다.
- 실제 메모리 회수 수단은 3.7절의 `replace_deleted`뿐입니다.

### 3.7 빈자리 재사용 `replace_deleted` (`L1037-L1091`)

생성자에서 `allow_replace_deleted=true`를 주고, `addPoint(..., replace_deleted=true)`로 호출해야 동작합니다.

```
lock(새 label 단위 잠금)
id = deleted_elements에서 아무 원소 하나를 꺼냄      # unordered_set::begin()
if 없음: 일반 삽입
else:
    old = label(id);  label(id) = new
    label_lookup_.erase(old);  label_lookup_[new] = id
    unmarkDeletedInternal(id)
    updatePoint(vec, id, 1.0)                        # 3.4절 + 3.5절
```

- 삭제된 노드의 슬롯, 레벨, 기존 간선을 그대로 재활용하고 `updatePoint`로 새 위치에 맞게 다시 연결합니다. 기존 간선이 남아 있으므로 3.4절에서 설명한 단방향 inbound 간선 문제가 그대로 생깁니다. 이전 벡터 근처의 노드들이 새 벡터를 가리키는 간선을 계속 가집니다.
- 코드 주석(`L1071`)에 "삭제된 원소에 대한 동시 작업은 없다고 가정한다"고 적혀 있습니다. 잠금은 **새** 라벨 기준으로만 잡으므로, 같은 순간에 **이전** 라벨로 `unmarkDelete`를 호출하면 경쟁 상태가 됩니다. 이 때문에 `unmarkDelete`는 replace 모드에서 안전하지 않다고 주석(`L994-L995`)에 명시되어 있습니다.
- 같은 label로 다시 `addPoint`를 호출하는데 그 원소가 삭제 상태이면, replace 모드에서는 오류를 반환합니다(`L1288-L1293`). replace 모드가 아니면 삭제 표시를 풀고 `updatePoint`를 실행합니다.

## 4. 동시성 모델

| 잠금 | 범위 | 보호 대상 |
|---|---|---|
| `label_op_locks_[label & 0xFFFF]` | 65,536개로 나눈 striped lock | 같은 label에 대한 삽입·삭제·갱신을 직렬화 |
| `label_lookup_lock` | 전역 1개 | `label → 내부 ID` 해시맵 |
| `link_list_locks_[id]` | 노드마다 1개 | 해당 노드의 모든 계층 인접 리스트 |
| `global` | 전역 1개 | `maxlevel_`과 진입점. 새 노드가 최고 레벨을 갱신할 때만 삽입이 끝날 때까지 유지 |
| `deleted_elements_lock` | 전역 1개 | 재사용 가능한 슬롯 집합 |

- 삽입 스레드는 자기 노드의 잠금을 삽입이 끝날 때까지 유지하고, 탐색하거나 역방향 간선을 추가할 때 다른 노드의 잠금을 잠깐씩 잡습니다. 그래서 삽입끼리는 노드 단위의 세밀한 병렬성을 가집니다.
- 검색 경로(`searchKnn`, `searchBaseLayerST`)는 **잠금을 전혀 잡지 않습니다.** README도 `add_items`와 `knn_query`는 서로 thread-safe하지 않다고 명시합니다. 온라인 서비스에서 읽기와 쓰기를 섞으려면 외부에서 reader-writer lock을 걸거나, 읽기용 인덱스와 쓰기용 인덱스를 따로 두고 교체하는 방식이 필요합니다.
- `resizeIndex`는 `realloc`으로 메모리를 다시 할당하므로 어떤 연산과도 동시에 실행할 수 없습니다.
- 용량 `max_elements`는 생성할 때 정합니다. 인덱스가 가득 차면 삽입이 실패하므로 resize나 replace가 필요합니다.

## 5. 실험

환경은 4코어 Xeon 2.8GHz, GPU 없음입니다. 데이터는 64차원 가우시안 혼합(군집 100개) 20,000점이고, `M=16`, `ef_construction=200`을 썼습니다. 정답은 살아 있는 점 전체에 대한 브루트포스 결과입니다. 원시 출력은 `experiments/hnsw_online_update/results.txt`, `delete_latency.txt`에 있습니다.

### 5.1 교체(churn)와 갱신 후의 recall@10

쿼리 500개, 4스레드로 측정했습니다. "새 인덱스"는 같은 시점의 살아 있는 점으로 처음부터 다시 만든 인덱스입니다.

| 시나리오 | ef=20 | ef=50 | ef=100 |
|---|---|---|---|
| 초기 구축 | 0.9744 | 0.9982 | 1.0000 |
| 30% 삭제 (tombstone) | 0.9914 | 0.9996 | 1.0000 |
| 30% 삭제 → 새 인덱스 | 0.9820 | 0.9994 | 1.0000 |
| 교체 5회 후 (누적 22,000점 교체) | 0.9588 | 0.9972 | 0.9992 |
| 교체 5회 → 새 인덱스 | 0.9730 | 0.9988 | 0.9998 |
| 이어서 20% 값 갱신 | 0.9600 | 0.9974 | 0.9994 |
| 20% 값 갱신 → 새 인덱스 | 0.9754 | 0.9988 | 0.9996 |

- tombstone만 있는 상태는 새 인덱스보다 오히려 recall이 높았습니다. 삭제된 노드가 경유지로 남아 있어 그래프가 더 촘촘하기 때문으로 보입니다.
- 교체와 값 갱신을 거친 인덱스는 `ef=20`에서 새 인덱스보다 약 1.5%p 낮았습니다. `ef≥50`에서는 차이가 0.002 이내였습니다. 실무에서는 ef를 조금 올리면 품질 저하를 상쇄할 수 있는 수준입니다.

삽입 처리량 (4스레드):

| 경로 | 처리량 | 점당 시간 |
|---|---|---|
| 신규 삽입 | 20,000점 / 0.98 s | 약 49 µs |
| `replace_deleted` 삽입 | 4,000점 / 1.06~1.15 s | 약 270 µs (약 5.5배) |
| 같은 label로 값 갱신 | 4,000점 / 1.04 s | 약 260 µs |

### 5.2 삭제 비율에 따른 검색 지연

쿼리 1,000개, 1스레드, `ef=20`으로 측정했으며, 5회 반복 측정의 중앙값입니다.

| 삭제 비율 | recall@10 | 쿼리당 지연 (µs) |
|---|---|---|
| 0% | 0.9746 | 31.0 |
| 10% | 0.9807 | 23.0 |
| 30% | 0.9904 | 27.9 |
| 50% | 0.9952 | 29.2 |
| 70% | 0.9978 | 39.4 |
| 90% | 0.9978 | 100.5 |

- 삭제 비율 50%까지는 지연이 거의 늘지 않았습니다. 70%부터 늘기 시작해 90%에서는 0%일 때의 약 3.2배가 되었습니다.
- 10~50% 구간의 지연이 0%보다 낮게 나온 것은 측정 잡음과 캐시 효과로 보입니다. 이 구간에서는 증가 추세가 뚜렷하지 않다고 해석하는 편이 안전합니다.
- recall이 오히려 올라가는 이유는 살아 있는 점이 줄어 문제가 쉬워지고, 같은 ef 안에서 더 많은 노드를 확장하기 때문입니다.
- 결론적으로 tombstone은 품질보다 **비용** 문제입니다. 삭제 비율이 높은 워크로드에서는 `replace_deleted`나 주기적인 재구축으로 tombstone 비율을 관리해야 합니다.

## 6. 연구 관점의 개선 여지

코드 분석과 실험에서 드러난 한계를 연구 주제 후보로 정리했습니다.

| 한계 | 위치 | 개선 방향 |
|---|---|---|
| 삭제된 노드가 계속 경유지로 쓰여 지연이 늘어남 | `searchBaseLayerST` | 삭제된 노드의 이웃을 우회 간선으로 이어 주는 국소 복구(FreshDiskANN의 delete consolidation 방식) |
| 단방향 inbound 간선이 갱신 후에도 남음 | `updatePoint` | 역방향 인접 리스트를 유지하거나, 주기적으로 오래된 간선을 정리 |
| 휴리스틱이 간선을 버려 inbound가 0인 노드가 생길 수 있음 | `mutuallyConnectNewElement` | `keepPrunedConnections` 구현 또는 inbound 하한 보장 |
| 삭제된 진입점을 계속 사용 | `addPointWithLevel`, `searchKnn` | 진입점이 삭제되면 살아 있는 최고 레벨 노드로 교체 |
| 갱신 비용이 삽입의 약 5.5배 | `updatePoint` | `updateNeighborProbability < 1`로 표본을 줄이고 recall 영향 측정 |
| 검색과 쓰기의 동시 실행 불가 | 전역 | 인접 리스트 버전 관리(RCU, copy-on-write)로 잠금 없는 읽기 |
| 물리 삭제·압축 없음 | 전역 | 내부 ID를 다시 매기는 오프라인 압축 루틴 |

## 부록: 실험 재현

```bash
git clone https://github.com/nmslib/hnswlib third_party/hnswlib   # 분석 커밋: ca42672
pip install numpy ./third_party/hnswlib
python3 experiments/hnsw_online_update/churn_recall.py
python3 experiments/hnsw_online_update/delete_latency.py
```
