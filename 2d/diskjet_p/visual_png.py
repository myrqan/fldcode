import numpy as np
import matplotlib.pyplot as plt
import os

# --- 設定パラメータ ---
params = {
    "dat_dir": "dat",         # Fortranデータがあるディレクトリ
    "png_dir": "png_output",  # PNGを出力するディレクトリ
    "var_name": "ro",         # 変数名
    "npx": 2, "npz": 2,       # ランク分割数
    "lx": 200, "lz": 300,     # 1ランクあたりの格子数
    "start_step": 0,          # 開始ステップ
    "end_step": 102,           # 終了ステップ
    "vmin": None,             # カラーバーの最小値 (Noneで自動)
    "vmax": None,             # カラーバーの最大値 (Noneで自動)
    "cmap": "jet"             # カラーマップ
}

def load_coordinates():
    """ x.dat, z.dat, t.dat を読み込む """
    try:
        x = np.fromfile('dat/x.dat', dtype=np.float64)
        z = np.fromfile('dat/z.dat', dtype=np.float64)
        t = np.fromfile('dat/t.dat', dtype=np.float64)
        return x, z, t
    except Exception as e:
        print(f"Error loading coordinates: {e}")
        return None, None, None

def load_aggregated_step(istep, p):
    """ 指定ステップの全ランクデータを統合する """
    global_data = np.zeros((p["npx"] * p["lx"], p["npz"] * p["lz"]), dtype=np.float64)
    
    for iz in range(p["npz"]):
        for ix in range(p["npx"]):
            fname = f"{istep:03d}_{p['var_name']}_{ix}_{iz}.dat"
            fpath = os.path.join(p["dat_dir"], fname)
            
            if os.path.exists(fpath):
                local_data = np.fromfile(fpath, dtype=np.float64)
                # Fortran orderでリシェイプ
                grid = local_data.reshape((p["lx"], p["lz"]), order='F')
                # グローバル配列に配置
                global_data[ix*p["lx"]:(ix+1)*p["lx"], iz*p["lz"]:(iz+1)*p["lz"]] = grid
                
    return global_data

def save_frames_as_png(p):
    # 1. 準備
    x, z, t = load_coordinates()
    if x is None: return

    # PNG出力ディレクトリの作成
    os.makedirs(p["png_dir"], exist_ok=True)
    print(f"Output directory: {p['png_dir']}")

    # 座標メッシュの作成
    X, Z = np.meshgrid(x, z, indexing='ij')

    # 2. ステップループ
    for istep in range(p["start_step"], p["end_step"] + 1):
        print(f"Processing step: {istep:03d}...")

        # データの読み込み
        data = load_aggregated_step(istep, p)
        
        # 時刻の取得
        current_time = t[istep] if istep < len(t) else t[-1]

        # 図の作成 (毎回新しく作る)
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # プロット (vmin, vmaxを指定しなければフレーム毎に自動調整される)
        pcm = ax.pcolormesh(X, Z, np.log10(data), cmap=p["cmap"], shading='auto',
                            vmin=p["vmin"], vmax=p["vmax"])
        
        # 装飾
        fig.colorbar(pcm, label=p["var_name"])
        ax.set_xlabel('Physical X [m]')
        ax.set_ylabel('Physical Z [m]')
        ax.set_title(f"Variable: {p['var_name']} | Step: {istep:03d} | Time: {current_time:.4f} s")
        ax.set_aspect('equal') # アスペクト比を物理座標に合わせる

        # PNG保存
        png_fname = os.path.join(p["png_dir"], f"frame_{istep:03d}.png")
        plt.savefig(png_fname, dpi=150, bbox_inches='tight')
        
        # メモリ開放 (非常に重要)
        plt.close(fig)

    print("All frames saved successfully.")

if __name__ == "__main__":
    # パラメータを調整して実行
    # 例: カラースケールを固定したい場合は vmin, vmax に数値を指定する
    # params["vmin"] = -1.0
    # params["vmax"] =  1.0
    
    save_frames_as_png(params)
