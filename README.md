# Droplets Simulation
Simulation of Virus Droplets Behavior in AFDET

## 使い方
  ※このブランチは **GNU Fortran, GNU Make** 用です。`master`ブランチへは絶対にマージしないでください。
  コンパイルに`make`コマンドを使います（`GNU make`のインストールが必要）。

  1. 「SampleCase」ディレクトリを複製したのち、名前を変更する（ケース名を付ける）。
  2. ケースディレクトリ内の条件ファイル(condition.nml, initial_position.csv)を編集。
  3. Makefileのあるディレクトリで `make` コマンド（コンパイル）。
  4. `.\bin\droplet.exe`で実行。ケース名を入力して計算開始。

## 条件ファイル(condition.nml, initial_position.csv)解説
  ### condition.nml
  - **リスタート位置 num_restart**
    - 通常は`0`を指定
    - `1以上`にすると、その値に対応するbackupファイルが読み込まれ、そこからリスタートが始まる
    - `-1`にすると、backupファイル(.bu)が読み込まれ、それを初期飛沫分布とする。backupファイル名は自由に指定可能。
  - **飛沫周期発生 preriodicGeneration**
    - 1秒当たりの発生飛沫数（整数）を指定
    - 初期配置飛沫をすべてNonActiveにしたのち、順次Activateしていくので、初期配置数が飛沫数の上限となる
  - **気流データファイル名 path2FlowFile**
    - 実行ディレクトリからの相対パス、もしくは絶対パスを指定
    - 現在可能な流れ場ファイル：
      - VTK
      - INP
      - FLD
    - CUBE格子(PLOT3D)は、予め非構造格子に変換してから計算してください。
    - .arrayファイルを指定する場合、別途メッシュファイルが必要なので、`meshFile = ***`と指定する
  - **ステップ数オフセット OFFSET**
    - 飛沫計算を、流体連番ファイルの途中の番号から始めたいときに指定
  - **気流データを周期的に用いる場合の先頭と末尾 LoopHead, LoopTail**
    - 任意の区間の流体連番ファイルを繰り返し用いるときに指定（例えば呼吸のサイクル）
    - `(先頭) = (末尾)` とすると、そのステップ数到達後は流れ場の更新が起こらなくなる
    - `(先頭) > (末尾)` とすれば、特殊な処理は起こらず、流体連番ファイルが順番に読み込まれる
  ### initial_position.csv
  - 初期飛沫の配置帯（直方体）を設定する
  - 左から順に、直方体の中心座標(x,y,z), 直方体の幅(x,y,z)
  - 改行すれば配置帯を複数設定できる

## 外部サブルーチン「dropletManagement」
  廃止しました。ボックス等で任意の場所の飛沫数をカウントしたい場合はdropletCount.f90を適宜書き換えて実行してください。

## 方程式

  ### 飛沫の蒸発方程式

  $$ \frac{dr}{dt} \space = \space -\left(1-\frac{RH}{100}\right) \cdot \frac{D e_{s}(T)}{\rho_{w} R_{v} T r} $$
  
  プログラム内では、２次精度ルンゲクッタ法で解いている。
  
  ### 飛沫の運動方程式

$$ m \frac{d \mathbf{v}}{dt} \space = \space m \mathbf{g} \space + \space C_D (\mathbf{v}) \space \cdot \space \frac{1}{2} \rho_a S \left | \mathbf{u}_a - \mathbf{v} \right | (\mathbf{u}_a - \mathbf{v}) $$

  プログラム内では、上式を無次元化・離散化した次式を解いている。
    
$$ \bar{\mathbf{v}}^{n + 1} \space = \space \frac{\bar{\mathbf{v}}^{n} \space + \space (\bar{\mathbf{g}} \space + \space C \bar{\mathbf{u}}_a)\Delta \bar{t}}{1 \space + \space C\Delta \bar{t}} \quad \left ( C \space = \space \frac{3 \rho_a}{8 \rho_w} \frac{C_D ( \mathbf{v}^{n} ) \left | \bar{\mathbf{u}}_a - \bar{\mathbf{v}}^{n} \right |}{\bar{r}^{n+1}} \right ) $$

## サブプログラム
  `make [subProgramName]`で実行ファイルを作成できる。
  - CUBE2USG
    - CUBE格子を、非構造格子に変換できる
  - droplet2CSV
    - 飛沫計算結果を再度読み込み、統計データ（浮遊数推移など）をCSVファイルに書き出す
  - dropletCount
    - 飛沫計算結果を再度読み込み、カウントボックスを通過した飛沫数を調べる。optionディレクトリ内の"boxList.csv"を、ケースディレクトリに配置する必要がある。
  - initialTranslate
    - 飛沫の初期配置データを読み込み、任意の座標への回転、平行移動を行う。**by Konishi**

## おまけ機能
  - **複数ケース連続実行**
    - 実行時にTXTファイル名を入力すると、そのファイルに列挙された複数ケースを連続実行できる
  - **basicSetting.nml**
    - optionディレクトリ内にある。付着判定のオンオフや、飛沫間合体の設定が可能。初期半径分布ファイルの指定も可能。

