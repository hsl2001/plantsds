# Multi-Pangenome SD 분석 프레임워크

## 1. 3×3 SD 진화 분류 체계 (Multi-Level Framework)

단일 유전체(Level 1)의 SD가 종 내(Level 2)와 종 간(Level 3)으로 확장될 때의 빈도와 보존성을 3×3 매트릭스로 맵핑하여 9가지 진화적 궤적을 정의합니다.

- **종 내 (Level 2) 기준**: Core (≥95%), Dispensable (5~95%), Private (≤5%)
- **종 간 (Level 3) 기준**: 7개 종에서의 보존 (Core), 일부 보존/다형성 (Dispensable), 단일 종 유일 (Private)

| 종 내(L2) \ 종 간(L3) | **Core (C) in Multi-Pan** | **Dispensable (D) in Multi-Pan** | **Private (P) in Multi-Pan** |
| :--- | :--- | :--- | :--- |
| **Core (C)**<br/>*(in Species)* | **C-C (초고대 필수 SD)**<br/>종 내 고정 & 식물계 전체 보존.<br/>*(기본 대사, 세포 분열 등)* | **C-D (계통 특이적 적응)**<br/>우리 종엔 고정, 타 종엔 다형성.<br/>*(환경에 완벽히 적응 완료)* | **C-P (종 특이적 혁신)**<br/>우리 종에만 유일하게 고정됨.<br/>*(새로운 형질/종분화의 핵심)* |
| **Dispensable (D)**<br/>*(in Species)* | **D-C (퇴화 중인 필수 SD)**<br/>타 종엔 다 고정인데 우리 종만 다형성.<br/>*(기능 소실 또는 구조적 퇴화)* | **D-D (수렴/평행 진화)**<br/>우리도 다형성, 타 종도 다형성.<br/>*(지속적 군비경쟁, NBS-LRR 등)* | **D-P (신흥 적응 SD)**<br/>우리 종에만 존재하며 확산 중.<br/>*(최근의 환경 적응 선택압)* |
| **Private (P)**<br/>*(in Species)* | **P-C (소실 직전의 흔적)**<br/>타 종엔 필수이나 우리 종엔 희귀.<br/>*(진화적 도태/멸종 직전)* | **P-D (일시적/불완전 SD)**<br/>타 종엔 다형성이나 우리는 희귀.<br/>*(유입/소실이 잦은 영역)* | **P-P (최근 발생 돌연변이)**<br/>우리 종 극소수 개체에만 갓 생겨남.<br/>*(De novo SD)* |

---

## 2. 데이터셋 출처 (7종, 총 750+ 게놈)

| # | 종 (Species) | 논문 출처 | 주요 다운로드 경로 및 Accession |
| :---: | :--- | :--- | :--- |
| 1 | **Tomato**<br/>(T2T 100+) | Shi et al. 2026, *Nat Genet* | - Zenodo: [17878268](https://zenodo.org/records/17878268)<br/>- CNCB: PRJCA030093<br/>- NCBI: PRJNA1201608 |
| 2 | **Grapevine**<br/>(70+) | Liu et al. 2024, *Nat Genet* | - Zenodo: [10851548](https://zenodo.org/records/10851548), [10846425](https://zenodo.org/records/10846425)<br/>- NCBI BioProject: PRJNA1018808 등 다수 |
| 3 | **Watermelon**<br/>(150+) | Sun et al. 2026, *Nat Genet* | - CuGenDBv2 (FTP): [graph_pangenome/assembly](http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/graph_pangenome/assembly/)<br/>- NCBI: PRJNA1272048 |
| 4 | **Citrus**<br/>(80+) | Huang et al. 2023, *Nat Genet* | - HZAU DB: [download.php](http://citrus.hzau.edu.cn/download.php)<br/>- Figshare: [a1e8071844912a7495ac](https://figshare.com/s/a1e8071844912a7495ac) |
| 5 | **Marchantia**<br/>(133+) | Beaulieu et al. 2025, *Nat Genet* | - MarpolBase: [pangenome_assemblies.tar.gz](https://marchantia.info/download/pangenome_assemblies.tar.gz)<br/>- NCBI: PRJNA931118, PRJNA1021402 |
| 6 | **Arabidopsis**<br/>(69) | 2024, *Nat Genet* | - Edmond Dataverse: [doi:10.17617/3.AEOJBL](https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.AEOJBL)<br/>- GitHub: [qclian/Pan_Ath](https://github.com/qclian/Pan_Ath) |
| 7 | **Rice**<br/>(149) | 2025, *Nature*<br/>([10.1038/s41586-025-08883-6](https://doi.org/10.1038/s41586-025-08883-6)) | - Figshare: [25697817](https://doi.org/10.25452/figshare.plus.25697817)<br/>- ENA/NCBI: PRJEB73710<br/>- NGDC: PRJCA024131<br/>- DB: RicePandb ([ricepandb.ncgr.ac.cn](http://ricepandb.ncgr.ac.cn)) |


---

## 3. Multi-Pangenome 3×3 활용 연구 계획

### 3.1 SD 검증 (Tool Verification)
- 실제 알려진 생물 종(Human 등)에서의 비율과 유사한 수준으로 SD가 감지되는지 검증
- 시뮬레이션 데이터를 활용하여 민감도및 정밀도 측정
- 특정 SD 서열을 전체 유전체 대상 BLAST 수행 시 본 툴에서 감지한 위치와 일치하는지 확인
- Eichler lab에서 개발한 assembly quality control pipeline 돌려보기

> **단기 목표**: 위 검증 단계를 완료하여 올해 추계 학회 발표 자료 준비

### 3.2 생물학적 타당성 검증
- $dN/dS$ 분석을 통한 선택압 확인:
  - `C-C` 그룹이 타 패턴 대비 음성 선택을 더 강하게 받는지 검증
  - `P-P` 그룹이 더 강한 양성 선택압을 받는지 확인
- 집단유전학적 지표 분석:
  - `P-P` 그룹에서 Selective sweep 패턴이 더 빈번하게 나타나는지 확인
  - $F_{ST}$, $\pi$ 등의 집단 유전 지표가 `P-P` 그룹에서 더 높게 나타나는지 검증

### 3.3 Functional Enrichment Analysis
- 9개 3×3 매트릭스 각 그룹에 속한 유전자 대상 Functional Enrichment Analysis 수행
- 서열 기반 분석: eggNOG-mapper를 활용한 기능 주석 분석
