# fldcode
流体計算コード
各モデルは`make`でコンパイル(gfortranを利用)，`./a.out`で実行，`python drawing.py`で描画ができる

## 1d/は1次元計算
- shocktube/は衝撃波管問題，これはうまくいってるはず．
- shocktube2/は衝撃波管問題を，Rubin&Burstin(1967)で示されている，式(7)あるいは式(10)（でソース項ゼロとする）に従って計算を行った．
- sedov/は1次元セドフ解を解いた．
- _sedovは．何かがうまく行かなかった時のやつ．ソース項の取り扱いだと思う．

## 2d/は2次元計算
- nan・prで負数が出ているので，それを解決すべき，
- xmax，zmaxのあたりからnanが増殖している気がする，

### to do
- 1d/sedovを頑張る
- $\partial u/\partial t + \partial F/\partial x = R(x, t)$のように右辺にsource termがある場合の取り扱いを考える．
    - modified Lax-Wendroff法の文献を調査する．→Rubin&Burstein(1967)に従って計算を行う．
    - CANSの計算方法と同じように，式(7)，式(10)を用いた計算で実行をしてみる．
下の式は$t + dt/2$での値を求める方法で，これは，CANSでは用いられていない．
```math
\frac{\partial u}{\partial t} + \frac{\partial F}{\partial x} = R(x, t)
```
```math
u_{j+1/2}^{n+1/2} = \frac{1}{2} \left[u_{j+1}^n + u_j^n \right] - \frac{\Delta t}{2 \Delta x} \left[F_{j+1}^n - F_j^n \right] + \frac{\Delta t}{2} \frac{1}{2} \left[ R_j^n + R_{j+1}^{n} \right]
```
```math
u_{j}^{n+1} = u_{j}^n - \frac{\Delta t}{\Delta x} \left[F_{j+1/2}^{n+1/2} - F_{j-1/2}^{n+1/2} \right] + \Delta t \frac{1}{2} \left[R_{j+1/2}^{n+1/2} + R_{j-1/2}^{n+1/2}\right]

```

- CFL条件を満たすようにdtの値をとるように（dtをstep毎に可変の値と）する

