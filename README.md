# droplets_simulation
Simulation of Droplets Behavior in AFDET

## grid_info.f90
  前処理用プログラム。格子隣接関係・境界面情報を出力。以下をInclude。
  - flow_field.f90
  - cases_reader.f90
  - csv_reader.f90
  
## virus_main.f90
  飛沫計算メインプログラム。以下をInclude。
  - flow_field.f90
  - virus_mod.f90
  - motion_mod.f90
  - cases_reader.f90
  - csv_reader.f90

## 使用するモジュール
  - flow_field.f90   :流れ場の格子データに関する変数・手続き集
  - virus_mod.f90    :飛沫に関する変数・手続き集
  - motion_mod.f90   :飛沫の運動に関する変数・手続き集
  - cases_reader.f90 :連続実行用ファイルcases.csvを読み込む
  - csv_reader.f90   :一般的なCSVファイルを読み込む

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m_d&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m_d&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\mathbf{v}^{n&plus;1}&space;=&space;\frac{\mathbf{v}^{n}&space;&plus;&space;(\mathbf{g}&space;&plus;&space;C\mathbf{u}_a)\Delta&space;t}{1&plus;C\Delta&space;t}" />
