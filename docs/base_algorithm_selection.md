# 온라인 갱신 가속 하드웨어 연구의 기준(base) 알고리즘 선정

근거: `reports/벡터 DB 온라인 갱신 알고리즘과 하드웨어.md`와 `research_notes/` 다섯 편, 그리고 `docs/hnswlib_online_update_analysis.md`의 코드 분석.

## 결론

**기준 알고리즘: Vamana(DiskANN) 계열의 in-place 갱신, 구체적으로 IP-DiskANN의 삽입·삭제 규칙을 OdinANN식 고정 크기 레코드 SSD 레이아웃 위에서 실행하는 구성.**

- 소프트웨어 기준 구현: `microsoft/DiskANN`의 dynamic index(FreshVamana + lazy delete + consolidate_deletes). 여기에 IP-DiskANN의 in-place 삭제(Algorithm 5·6)를 얹는다. IP-DiskANN 원 구현은 비공개 Rust이므로 직접 구현해야 하지만, 알고리즘은 삽입 커널의 재사용이라 구현량이 작다.
- 저장 모델: OdinANN(FAST '26)의 "벡터 + 이웃 ID R개"를 담은 고정 크기 레코드, 페이지당 여유 슬롯으로 GC 없이 갱신을 묶는 방식.
- 평가 프로토콜: NeurIPS BigANN 2023 streaming track의 runbook(SlidingWindow, ExpirationTime, Clustered)과 MSTuring/Wiki-Cohere 데이터셋. 비교 대상은 같은 계열인 FreshDiskANN, IP-DiskANN, OdinANN, Yi, NAVIS, 그리고 GPU 쪽 GrAND/SVFusion.

보조 기준(두 번째 후보): **SPFresh/LIRE**(SOSP '23, 오픈소스). 연구 대상이 SmartSSD처럼 "장치 내부 연산"이라면 posting = 블록 묶음이라는 대응이 자연스러워 두 번째 기준으로 둘 가치가 있다. 다만 아래 이유로 1순위는 아니다.

## 선정 이유

### 1. 갱신의 핵심 커널이 검색과 같다

IP-DiskANN 계열에서 삽입과 삭제는 모두 다음 세 단계로 이루어진다.

| 단계 | 삽입 | 삭제 (in-place) | 비용 성격 |
|---|---|---|---|
| ① position seeking | 새 벡터로 GreedySearch(L≈75~128) | 삭제 벡터로 GreedySearch(l_d=128, k=50) | 무작위 4 KB 읽기 약 100회 + 거리 계산 약 8천 회. OdinANN 기준 갱신 지연의 최대 85% (NAVIS 측정) |
| ② 이웃 선택 | RobustPrune(α=1.2, R) | 방문 집합 중 삭제 노드를 가리키는 노드를 in-neighbor로 근사하고, 각 노드에 대체 간선 c=3개 선택 | 후보 집합 내부의 거리 계산 O(L·R), 규칙적이고 작은 연산 |
| ③ 간선 패치 | 이웃 최대 R개의 레코드에 역방향 간선 추가, 차수 초과 시 prune | in-neighbor 레코드 수정, 주기적으로 거리 계산 없는 dangling-edge 정리 | 레코드 단위 read-modify-write, 쓰기 증폭의 원인 |

①은 검색 그 자체다. NDSearch(ISCA '24), CXL-ANNS(ATC '23), Cosmos(CAL '25)처럼 이미 검증된 검색용 데이터패스를 그대로 가져와 갱신에 쓸 수 있다. 하드웨어 연구자 입장에서는 "검색 가속기를 갱신까지 확장"이라는 명확한 논제가 생기고, ②와 ③이 새로 설계할 부분이 된다.

SPFresh/LIRE는 반대다. 핫패스(posting 끝에 append)는 이미 1.5 ms로 싸고, 비싼 분할·재배치는 삽입의 0.4%에서만 일어난다. 희귀한 경로를 가속해서는 전체 이득이 작다.

### 2. 재현율이 안정적이라 하드웨어 효과를 깨끗하게 측정할 수 있다

- α>1 재가지치기 덕에 50회 이상의 5~50% 교체 주기에서도 재현율이 유지된다(FreshDiskANN Fig. 2, IP-DiskANN Table 1).
- in-place 설계라 병합 중 지연 급등(DiskANN P99.9 > 20 ms, OdinANN 측정 1.54~2.44배 변동)이 없다. 가속기의 효과가 병합 스파이크에 묻히지 않는다.
- SPFresh는 군집 단위가 거칠어 시작 재현율이 0.68~0.75라는 독립 보고(LSM-VEC, Yi)가 있다. 재현율 상한이 낮으면 "가속해도 품질은 그대로"라는 주장을 하기 어렵다.

### 3. 비교 기준선과 벤치마크가 같은 계열 안에 갖춰져 있다

FreshDiskANN(2021) → IP-DiskANN(2025) → OdinANN(FAST '26) → Yi/NAVIS/MERIT(2026) → GrAND(GPU, 2026)이 모두 Vamana 기반이고 BigANN streaming runbook으로 평가한다. 같은 runbook을 돌리면 소프트웨어 최신 결과와 직접 비교할 수 있다. Azure Cosmos DB가 같은 in-place 삭제를 채택했으므로 산업적 타당성도 있다.

### 4. 하드웨어 공백이 정확히 이 커널에 있다

조사한 어떤 출처도 SSD·CXL·PIM 장치 내부에서 ①~③을 수행하지 않는다(보고서 5~6절). 소프트웨어 쪽이 io_uring, SPDK, 레이아웃 분리까지 다 써 봤는데도 position seeking이 호스트 CPU에 남아 있다는 것이 OdinANN·NAVIS·Yi가 공통으로 지목한 병목이다.

## 제외한 후보와 이유

| 후보 | 제외 이유 |
|---|---|
| hnswlib(HNSW) | 메모리 상주 전용. 삭제는 톰스톤뿐이고 replace 경로는 삭제당 O(R³) 간선을 만든다(IP-DiskANN 분석). in-neighbor 복구가 없어 3,000회 교체 후 도달 불가 노드 3~4%(MN-RU). 계층 구조가 장치 레이아웃을 복잡하게 한다. 메모리 내 참조 구현으로만 유지한다. |
| FreshDiskANN의 batch consolidation | 삭제당 O(R²) 간선과 전체 인덱스 스캔. 병합 사이에 재현율이 내려간다. 후속 연구가 이미 대체했다. |
| SPFresh/LIRE | 위 1·2항. 다만 in-storage 연구라면 보조 기준으로 유지. |
| GPU 동적 인덱스(Jasper, GrAND, SVFusion) | 2026년 프리프린트 중심, HBM 용량에 묶여 GPU당 1~2억 벡터, 배치 단위 가시성. 경쟁 상대로 두되 기준으로 삼지 않는다. |
| LSM-VEC, Ada-IVF, UBIS | 갱신 처리량 개선은 있으나 1M~100M 단일 노드 평가이고 커널이 ①~③과 다르다. |

## 하드웨어 설계로 이어지는 분해

| 커널 | 가속 위치 후보 | 참고할 기존 설계 |
|---|---|---|
| ① position seeking (무작위 읽기 + 거리 계산) | SSD 컨트롤러/FPGA 또는 CXL 장치 내 탐색 엔진 | NDSearch의 search-page 엔진, CXL-ANNS·Cosmos의 거리 오프로드 |
| ② RobustPrune / 대체 간선 선택 (후보 집합 내부 all-pairs 거리) | 장치 내 소형 거리 배열(systolic) | ANNA·ANSMET의 거리 유닛 |
| ③ 역방향 간선 패치 (레코드 RMW) | FTL 수준의 페이지 내 패치, OdinANN의 update combining을 장치에서 수행 | OdinANN 고정 레코드 + 여유 슬롯 |
| 삭제 후 dangling-edge 정리 (거리 계산 없는 순차 스캔) | in-storage 순차 스캔 | SmartSSD 스캔 커널 |

호스트는 신선 계층(in-memory 버퍼)과 스케줄링만 맡고, 장치가 ①~③을 처리하는 분할이 가장 자연스럽다.

## 첫 실험 제안

1. `microsoft/DiskANN` dynamic index에 IP-DiskANN 삭제를 구현하고 BigANN streaming runbook(MSTuring-10M/30M)을 돌린다.
2. 갱신 한 건당 ①~③의 시간·I/O·거리 계산 횟수를 분해 측정해 NAVIS의 85% 주장을 자체 환경에서 재확인한다.
3. 이 분해가 하드웨어 설계의 Amdahl 상한을 정한다. ①이 80% 이상이면 검색 엔진 재사용만으로도 큰 이득이 있고, ③이 크면 페이지 패치 엔진이 필요하다.
