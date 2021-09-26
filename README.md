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
  - drop_motion.f90
  - equation_mod.f90
  - cases_reader.f90
  - csv_reader.f90

## 使用するモジュール
  - flow_field.f90   :流れ場の格子データに関する変数・手続き集
  - drop_motion.f90  :飛沫の挙動に関する変数・手続き集
  - equation_mod.f90 :方程式系を扱うモジュール
  - cases_reader.f90 :連続実行用ファイルcases.csvを読み込む
  - csv_reader.f90   :一般的なCSVファイルを読み込む

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を無次元化・離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\bar{\mathbf{v}}^{n&plus;1}&space;=&space;\frac{\bar{\mathbf{v}}^{n}&space;&plus;&space;(\bar{\mathbf{g}}&space;&plus;&space;C\bar{\mathbf{u}}_a)\Delta&space;\bar{t}}{1&plus;C\Delta&space;\bar{t}}" />
