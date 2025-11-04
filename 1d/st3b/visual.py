import numpy as np
import matplotlib.pyplot as plt
import read

datdir = 'dat/'

# --- データ読み込み ---
tt = read.rd_0d(datdir+'t.dat')
xx = read.rd_0d(datdir+'x.dat')
ro = read.rd_1d(datdir+'ro.dat')
vx = read.rd_1d(datdir+'vx.dat')
vy = read.rd_1d(datdir+'vy.dat')
bx = read.rd_1d(datdir+'bx.dat')
by = read.rd_1d(datdir+'by.dat')
pr = read.rd_1d(datdir+'pr.dat')
# ei = read.rd_1d(datdir+'ei.dat')

# --- 描画設定 ---
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'

frames = [0, 5, 10, 14]   # 出力する時刻インデックス
frames = [i for i in range(15)]

for n in frames:
    time = tt[n]
    fig, axs = plt.subplots(3, 2, figsize=(8, 9), sharex=True)
    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    fig.suptitle(r'$t = $' + str(time)[:5], fontsize=14)

    # 軸を1次元配列に（操作しやすく）
    ax = axs.ravel()

    # --- 各パネル描画 ---
    ax[0].plot(xx, ro[n], 'k')
    ax[0].set_ylabel(r'$\rho$ (density)')
    ax[0].set_ylim(0.0, 1.5)

    ax[1].plot(xx, pr[n], 'r')
    ax[1].set_ylabel(r'$p$ (pressure)')
    ax[1].set_ylim(0.0, 1.5)

    ax[2].plot(xx, vx[n], 'b')
    ax[2].set_ylabel(r'$v_x$')
    ax[2].set_ylim(-0.5, 1.0)

    ax[3].plot(xx, vy[n], 'g')
    ax[3].set_ylabel(r'$v_y$')
    ax[3].set_ylim(-2.0, 0.5)

    ax[4].plot(xx, by[n]/np.sqrt(4*np.pi), 'm')
    ax[4].set_ylabel(r'$B_y/\sqrt{4\pi}$')
    ax[4].set_ylim(-1.5, 1.5)
    ax[4].set_xlabel(r'$x$')

    # 6番目のサブプロット（右下）は空白にしておくか、凡例などを表示
    ax[5].plot(xx, by[n]/np.sqrt(4*np.pi))
    ax[5].plot(xx, ro[n])
    ax[5].plot(xx, pr[n])
    ax[5].plot(xx, vx[n])
    ax[5].plot(xx, vy[n])

    ax[5].set_ylim(-2.0, 2.0)
    ax[5].set_xlabel(r'$x$')
    #ax[5].axis('off')
    #ax[5].text(0.5, 0.5, 'Brio & Wu (1988)\n1D MHD shock tube',
    #           ha='center', va='center', fontsize=12)

    # 全体設定
    for a in ax[:5]:
        a.set_xlim(0, 1)

    savname = f'graph/{str(n).zfill(3)}.png'
    plt.savefig(savname, dpi=300, bbox_inches='tight')
    plt.close(fig)

