# fldcode
流体計算コード
各モデルは`make`でコンパイル(gfortranを利用)，`./a.out`で実行，`python drawing.py`で描画ができる

## 1d/は1次元計算
- shocktube/は衝撃波管問題

### to do
- 1d/sedovを頑張る
- df/dt + dU/dx = S(x, t)のように右辺にsource termがある場合の取り扱いを考える．
    - modified Lax-Wendroff法の文献を調査する．
- CFL条件を満たすようにdtの値をとるように（dtをstep毎に可変の値と）する．
- 1d/sedov/main.f90 L 136 改修中

