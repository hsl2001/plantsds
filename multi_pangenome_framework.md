# Multi-Pangenome SD 분석 프레임워크

## 1. 3×3 SD 진화 분류 체계 (Multi-Level Framework)

단일 유전체(Level 1)의 SD가 **종 내(Level 2)**와 **종 간(Level 3)**으로 확장될 때의 빈도와 보존성을 3×3 매트릭스로 맵핑하여 9가지 진화적 궤적을 정의합니다.

- **종 내 (Level 2) 기준**: Core(≥95%), Dispensable(5~95%), Private(≤5%)
- **종 간 (Level 3) 기준**: 8개 종에서의 보존(Core), 일부 보존/다형성(Dispensable), 단일 종 유일(Private)

| 종 내(L2) \ 종 간(L3) | **Core (C) in Multi-Pan** | **Dispensable (D) in Multi-Pan** | **Private (P) in Multi-Pan** |
|:---|:---|:---|:---|
| **Core (C)**<br/>*(in Species)* | **C-C (초고대 필수 SD)**<br/>종 내 고정 & 식물계 전체 보존.<br/>*(기본 대사, 세포 분열 등)* | **C-D (계통 특이적 적응)**<br/>우리 종엔 고정, 타 종엔 다형성.<br/>*(환경에 완벽히 적응 완료)* | **C-P (종 특이적 혁신)**<br/>우리 종에만 유일하게 고정됨.<br/>*(새로운 형질/종분화의 핵심)* |
| **Dispensable (D)**<br/>*(in Species)* | **D-C (퇴화 중인 필수 SD)**<br/>타 종엔 다 고정인데 우리 종만 다형성.<br/>*(기능 소실 또는 구조적 퇴화)* | **D-D (수렴/평행 진화)**<br/>우리도 다형성, 타 종도 다형성.<br/>*(지속적 군비경쟁, NBS-LRR 등)* | **D-P (신흥 적응 SD)**<br/>우리 종에만 존재하며 확산 중.<br/>*(최근의 환경 적응 선택압)* |
| **Private (P)**<br/>*(in Species)* | **P-C (소실 직전의 흔적)**<br/>타 종엔 필수이나 우리 종엔 희귀.<br/>*(진화적 도태/멸종 직전)* | **P-D (일시적/불완전 SD)**<br/>타 종엔 다형성이나 우리는 희귀.<br/>*(유입/소실이 잦은 영역)* | **P-P (최근 발생 돌연변이)**<br/>우리 종 극소수 개체에만 갓 생겨남.<br/>*(De novo SD)* |

---

## 2. 데이터셋 출처 (8종, 총 768+ 게놈)

| # | 종 (Species) | 논문 출처 | 주요 다운로드 경로 및 Accession |
|:---|:---|:---|:---|
| 1 | **Tomato**<br/>(T2T 100+) | Shi et al. 2026, *Nat Genet* | - Zenodo (직접 다운로드): `https://zenodo.org/records/17878268`<br/>- CNCB: PRJCA030093, NCBI: PRJNA1201608 |
| 2 | **Grapevine**<br/>(70+) | Liu et al. 2024, *Nat Genet* | - Zenodo (직접 다운로드): `https://zenodo.org/records/10851547`, `10846425`<br/>- NCBI BioProject: PRJNA1018808 등 다수 |
| 3 | **Watermelon**<br/>(150+) | Sun et al. 2026, *Nat Genet* | - CuGenDBv2 (FTP): `http://cucurbitgenomics.org/v2/ftp/pan-genome/watermelon/`<br/>- NCBI: PRJNA1272048 |
| 4 | **Citrus**<br/>(80+) | Huang et al. 2023, *Nat Genet* | - HZAU DB (직접 다운로드): `http://citrus.hzau.edu.cn/download.php`<br/>- Figshare: `https://figshare.com/s/a1e8071844912a7495ac` |
| 5 | **Marchantia**<br/>(133+) | Beaulieu et al. 2025, *Nat Genet* | - MarpolBase (직접 다운로드): `https://marchantia.info/download/`<br/>- NCBI: PRJNA931118, PRJNA1021402 |
| 6 | **Cucumber**<br/>(100+) | Guan et al. 2026, *Nat Genet* | - NGDC GSA: PRJCA038097, PRJCA043228, PRJCA038675 |
| 7 | **Arabidopsis**<br/>(69) | 2024, *Nat Genet* | - GitHub Repository: `https://github.com/qclian/Pan_Ath` |
| 8 | **Rice**<br/>(149) | 2025, *Nature*<br/>(10.1038/s41586-025-08883-6) | - Figshare (직접 다운로드): `https://doi.org/10.25452/figshare.plus.25697817`<br/>- ENA/NCBI: PRJEB73710, NGDC: PRJCA024131<br/>- DB: RicePandb (`http://ricepandb.ncgr.ac.cn`) |
