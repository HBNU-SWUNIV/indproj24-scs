# 한밭대학교 SW중심대학 산학연계프로젝트 - 위성통신의 성능지표 시각화 시뮬레이션 연구

## **팀 구성**
### 지도교수
 - 방인규 교수님

### 기업체 
 - 이석근 대표

### 참여학생
 - 20221095 민현선 
 - 20221110 이다은

## Project Background
- ### 필요성
  - 넓은 범위와 지리적 제약을 극복할 수 있는 6G 위성통신 기술 필요
  - 저궤도 위성의 낮은 지연 시간과 넓은 커버리지 특성을 활용한 실시간 SNR 분석 필요
 
- ### 기존 해결책의 문제점
  - 통신 음영 지역 및 재난 상황에서 지상 5G 통신망의 한계 발생
  - 군집 모델, 궤도 평면의 수, 평면당 위성 수, 등 여러가지 파라미터를 조정하여 SNR 커버리지 맵 최적화 필요
  
## System Design
  - ### System Requirements
    <img width="678" alt="image" src="https://github.com/user-attachments/assets/372dcd20-966a-4430-8804-6fac2fa7c32f">
    1. 시간 및 위성의 위치 정보를 바탕으로 각 지역의 수신 SNR 값을 계산
    2. 수신 장치에서 GADM(Global Administrative Areas)에서 제공하는 대한민국의 지리 좌표 데이터를 기반으로 위도와 경도를 계산하여 각 좌표에 해당하는 SNR 값 매칭
    3. MATLAB 환경에서 데이터 처리 속도를 높이기 위해 샘플링 후 데이터 분석 수행

  - ### Scenario
    1) 실험① 에서 저궤도 위성 군집 모델 선정 (Walker Star 모델,  Waaler Delta 모델)
    2) 선정된 군집모델을 통해 실험② 에서 위성 배치별 SNR 비교
    3) 앞서 실험한 결과를 통해 적합한 군집 모델과 위성 배치 설정
    4) 설정한 조건들로 특정 위치의 SNR 커버리지맵 추출
    5) 특정 위치의 SNR(dBm)값 확인
    
## Case Study
본 연구는 저궤도 위성을 기반으로 대한민국 지역에 대한 SNR 커버리지 맵을 시뮬레이션하고 분석하기 위해 연구를 진행하였다. 
  - ### 저궤도 위성 모델 비교
    저궤도 위성 모델인 Walker Star와 Walker Delta 모델을 비교하였다.
    그 결과 **Walker Star 모델**을 사용한 경우 더 균등한 SNR Coverage Map이 생성되었다.
  - ### 위성 배치 별 SNR 커버리지 맵 비교
    위성의 배치를 균등한 경우와 균등하지 않은 경우를 비교하였다.
    그 결과 **위성의 배치가 균등하고 위성의 개수가 많을수**록 SNR Coverage Map 신호의 강도가 균일하게 나타났다.

=> Walker Star 모델을 사용하고 위성의 배치를 균등하게 배치하고, 위성의 개수가 많을 수록 더 좋은 SNR 커버리지 맵 출력이 가능하다. 
  
## Conclusion
Walker Star 모델과 균등한 위성 배치는 신호 강도와 통신 품질을 효과적으로 향상시킬 수 있음을 확인하였다. 
  -  균등한 위성 배치와 모델 선택은 SNR 커버리지 맵의 품질을 결정하는 주요 요소이다.
  -  기존 landarea.shp 파일에 존재하지 않는 국가의 SNR 커버리지 맵 출력이 가능하다.
  -  사용자가 원하는 지역(대전)의 SNR 추출이 가능하다.
  
## Project Outcome
- ### 2024년 추계종합학술발표회  
