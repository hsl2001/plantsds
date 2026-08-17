# Duplication map 관련 아이디어

### 코드에 추가해야 할 사항
- `-c` 옵션을 추가해서 "한 유전체(한 파일) 내부의 copy가 몇 개 이상이어야 레포트할지"를 넣어야 polyploid 대할 때 쓸 수 있음. 
- `-c` 가 2가 디폴트, 그러니까 한 유전체 내부에서 segmental duplication을 찾는 툴
- `-c` 가 1이면 duplication map을 만들 때 활용 가능. 
- `-c` 가 3이면 한 유전체 안의 카피가 최소 3개여야 레포트. 
- subclustering은 필요 없음. 
- 메모리 추가 절약해야 함.
- 속도가 지금보다 훨씬 빨라야 함.  
- 정확도 희생하면 안 됨. 
- 코드량이 더 늘어나면 안 됨. (subclustering 없어진 거 감안하더라도)


## 아이디어
이 방법(segtrace)을 이용하면
1) Pangenome에서 segmental duplication을 찾는 기법을 만들 수 있겠다. 
2) 식물 유전체 전체 (적어도 angiosperm)에서 duplication map을 그려볼 수 있겠다. 

### Duplication map
- Duplication map의 정의
    - 특정 윈도우로 타깃 유전체들을 전부 쪼개고, 
    - 그 윈도우들간의 all-versus-all 비교로 유사도를 분석한 다음
    - 그것으로 그래프를 그리는 것. 
    - 노드는 각 유전체의 각 윈도우, 엣지는 유사도 자체 (또는 특정 유사도 이상일 경우 1과 0)
    - 실로 거대한 크기의 그래프. 

Useful resources, technical advances, solving mechanisms 각각의 측면에서 강력한 당위성이 있는 어떠한 질문
> Duplication map이라는 개념이 다른 사람들이 이걸로 논문 안 쓰고는 못 배길 개념이 되도록.
- 무릇 좋은 개념이란, 그 개념에 대한 조사만으로도 어떤 논문이 성립 가능해야 함. 

plant-plant horizontal transfer를 연구하기 위한 자료가 될 수 있음. => 이게 뭐가 그렇게 중요함?


#### Useful resources
1) 당연하게도 angiosperm 한 수백개 되는 거의 duplication map을 만들 수 있음

#### Technical advances
1) Pangenome에서 SD를 찾는 기법
2) 다량의 유전체에서 duplication map을 그릴 수 있는 기법 => 무슨 쓸모?
- duplication map의 쓸모를 찾아내면 technical advance로 인정이 됨. 

#### Solving mechaism
- 아직 전혀 모르겠음. 