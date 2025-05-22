# fldcode
流体計算コード
各モデルは`make`でコンパイル(gfortranを利用)，`./a.out`で実行，`python drawing.py`で描画ができる

## 1d/は1次元計算
- shocktube/は衝撃波管問題，これはうまくいってるはず．

### to do
- 1d/sedovを頑張る
- $\partial u/\partial t + \partial F/\partial x = R(x, t)$のように右辺にsource termがある場合の取り扱いを考える．
    - modified Lax-Wendroff法の文献を調査する．
```math
\frac{\partial u}{\partial t} + \frac{\partial F}{\partial x} = R(x, t)\\
u_{j+1/2}^{n+1/2} = \frac{1}{2} \left[u_{j+1}^n + u_j^n \right] - \frac{\Delta t}{2 \Delta x} \left[F_{j+1}^n - F_j^n \right] + \frac{\Delta t}{2} \frac{1}{2} \left[ R_j^n + R_{j+1}^{n} \right]\\
u_{j}^{n+1} = u_{j}^n - \frac{\Delta t}{\Delta x} \left[F_{j+1/2}^{n+1/2} - F_{j-1/2}^{n+1/2} \right] + \Delta t \frac{1}{2} \left[R_{j+1/2}^{n+1/2} + R_{j-1/2}^{n+1/2}\right]

```

- CFL条件を満たすようにdtの値をとるように（dtをstep毎に可変の値と）する

