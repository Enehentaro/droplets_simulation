# droplets_simulation
Simulation of Droplets Behavior in AFDET

## main.f90
  飛沫計算メインプログラム。前処理プログラムは取り込まれました。
  実行時にファイルの有無から判断して、必要であれば前処理が勝手に行われます。
  現在可能な流れ場ファイル：
  - VTK
  - INP
  - FLD
  - PLOT3D

## 使い方
  コンパイルに`make`コマンドを使います（makeのインストールが必要）。
  - main.f90の冒頭部で、OSを指定する箇所があるので、適宜編集
  - Makefileを開き、コンパイラおよびコンパイルオプションを適宜編集（変数名：FC, FCFLAG）
  - ソースファイルのあるディレクトリで `make` コマンド実行
  - 「droplet」という実行ファイルが現れるので、それを実行

## 方程式

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を無次元化・離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\bar{\mathbf{v}}^{n&plus;1}&space;=&space;\frac{\bar{\mathbf{v}}^{n}&space;&plus;&space;(\bar{\mathbf{g}}&space;&plus;&space;C\bar{\mathbf{u}}_a)\Delta&space;\bar{t}}{1&plus;C\Delta&space;\bar{t}}" />
