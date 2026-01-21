import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

# カラーマップの準備
cmap_obj = sns.color_palette("coolwarm", as_cmap=True)

# --- 設定パラメータ ---
# 物理量や描画範囲に関する設定をここに集約しています
params = {
    "dat_dir": "dat",
    "png_dir": "png_output",
    "npx": 2, "npz": 2,
    "lx": 250, "lz": 250,
    "start_step": 0,
    "end_step": 15,
    
    # 背景（pcolormesh）の設定
    "base_var": "ro",
    "vmin": -4.0,
    "vmax": 1.0,
    "cmap": cmap_obj,
    
    # 速度ベクトル（quiver）の設定
    "plot_velocity": True,
    "vel_vars": ["vx", "vz"],
    "vel_skip": 50,        # 矢印の間隔（格子点数）
    "vel_scale": 1.0,      # 矢印の長さスケール
    "vel_color": "white",
    
    # 等高線（contour）の設定
    "plot_contour": True,
    "contour_var": "ay",
    "c_start": 0.0,
    "c_end": 0.5,
    "c_step": 0.01,
    "c_color": "white",
    "c_linewidth": 0.6
}

def load_coordinates():
    """ 座標軸データの読み込みを行います """
    try:
        x = np.fromfile('dat/x.dat', dtype=np.float64)
        z = np.fromfile('dat/z.dat', dtype=np.float64)
        t = np.fromfile('dat/t.dat', dtype=np.float64)
        return x, z, t
    except Exception as e:
        print(f"Error loading coordinates: {e}")
        return None, None, None

def load_aggregated_step(istep, var_name, p):
    """ 指定された変数の全ランクデータを統合して一つの二次元配列を返します """
    global_data = np.zeros((p["npx"] * p["lx"], p["npz"] * p["lz"]), dtype=np.float64)
    
    for iz in range(p["npz"]):
        for ix in range(p["npx"]):
            fname = f"{istep:03d}_{var_name}_{ix}_{iz}.dat"
            fpath = os.path.join(p["dat_dir"], fname)
            
            if os.path.exists(fpath):
                local_data = np.fromfile(fpath, dtype=np.float64)
                # Fortran順序でリシェイプし、グローバル配列の該当箇所に代入
                grid = local_data.reshape((p["lx"], p["lz"]), order='F')
                global_data[ix*p["lx"]:(ix+1)*p["lx"], iz*p["lz"]:(iz+1)*p["lz"]] = grid
    return global_data

def save_frames_as_png(p):
    """ メインの描画ループ処理です """
    x, z, t = load_coordinates()
    if x is None: return

    os.makedirs(p["png_dir"], exist_ok=True)
    # 座標メッシュの作成（indexing='ij'によりFortran配列と整合）
    X, Z = np.meshgrid(x, z, indexing='ij')

    # 各ステップの処理
    for istep in range(p["start_step"], p["end_step"] + 1):
        print(f"Processing step: {istep:03d}...")

        # 1. 背景用データの読み込みと描画
        base_data = load_aggregated_step(istep, p["base_var"], p)
        current_time = t[istep] if istep < len(t) else t[-1]

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.rcParams['font.size'] = 11
        plt.rcParams['font.family'] = 'Avenir'
        plt.rcParams['mathtext.fontset'] = 'stix'
        
        # 背景のカラーマップ描画（対数スケールを適用）
        pcm = ax.pcolormesh(
            X, Z, np.log10(base_data),
            cmap=p["cmap"],
            vmin=p["vmin"],
            vmax=p["vmax"],
            shading='auto'
        )
        
        # 2. 速度ベクトルのオーバーレイ
        if p["plot_velocity"]:
            vx = load_aggregated_step(istep, p["vel_vars"][0], p)
            vz = load_aggregated_step(istep, p["vel_vars"][1], p)
            
            # 格子点が多い場合は間引いて描画
            s = p["vel_skip"]
            ax.quiver(
                X[::s, ::s], Z[::s, ::s], 
                vx[::s, ::s], vz[::s, ::s],
                color=p["vel_color"],
                alpha=0.6,
                scale=p["vel_scale"],
                pivot='mid'
            )

        # 3. 等高線のオーバーレイ（ご指定の間隔設定を適用）
        if p["plot_contour"]:
            ay_data = load_aggregated_step(istep, p["contour_var"], p)
            # 指定された範囲と間隔でレベルを生成
            levels = np.arange(p["c_start"], p["c_end"], p["c_step"])
            
            ax.contour(
                X, Z, ay_data,
                levels=levels,
                colors=p["c_color"],
                linewidths=p["c_linewidth"],
                alpha=0.8
            )

        # グラフの装飾
        fig.colorbar(pcm, label=f"log10({p['base_var']})")
        ax.set_xlabel('X')
        ax.set_ylabel('Z')
        ax.set_title(f"Step: {istep:03d} | Time: {current_time:.4f} s")
        ax.set_aspect('equal')

        # 保存とメモリ解放
        png_fname = os.path.join(p["png_dir"], f"frame_{istep:03d}.png")
        plt.savefig(png_fname, dpi=168, bbox_inches='tight')
        plt.close(fig)

    print("Successfully completed.")

if __name__ == "__main__":
    save_frames_as_png(params)
