# droplets_simulation
Simulation of Droplets Behavior in AFDET

## grid_info.f90
  前処理用プログラム。
  
## virus_main.f90
  飛沫計算メインプログラム。以下をInclude。
  - flow_field.f90
  - virus_mod.f90
  - motion_mod.f90

  解くべき方程式は次の通り。  
<img src="https://latex.codecogs.com/gif.latex?m_d&space;\frac{d&space;\mathbf{v}}{dt}&space;=&space;m_d&space;\mathbf{g}&space;&plus;&space;C_D&space;\cdot&space;\frac{1}{2}\rho_a&space;S&space;\left&space;|&space;\mathbf{u}_a&space;-&space;\mathbf{v}&space;\right&space;|(\mathbf{u}_a&space;-&space;\mathbf{v})" />

  プログラム内では、上式を離散化した次式を解いている。  
<img src="https://latex.codecogs.com/gif.latex?\mathbf{v}^{n&plus;1}&space;=&space;\frac{\mathbf{v}^{n}&space;&plus;&space;(\mathbf{g}&space;&plus;&space;C\mathbf{u}_a)\Delta&space;t}{1&plus;C\Delta&space;t}" />
