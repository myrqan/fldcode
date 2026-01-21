import numpy as np
import matplotlib.pyplot as plt
import os

def load_coordinates():
    """
    x.dat, z.dat, t.dat から物理座標と時刻データを読み込む
    """
    try:
        # Fortranの real(8) は np.float64
        x_coords = np.fromfile('dat/x.dat', dtype=np.float64)
        z_coords = np.fromfile('dat/z.dat', dtype=np.float64)
        times    = np.fromfile('dat/t.dat', dtype=np.float64)
        return x_coords, z_coords, times
    except FileNotFoundError as e:
        print(f"Error: Coordinate file not found. {e}")
        return None, None, None

def load_aggregated_step(istep, var_name, npx, npz, lx, lz, dat_dir="dat"):
    """
    指定ステップの全ランクデータを統合する
    """
    global_data = np.zeros((npx * lx, npz * lz), dtype=np.float64)
    
    for iz in range(npz):
        for ix in range(npx):
            fname = f"{istep:03d}_{var_name}_{ix}_{iz}.dat"
            fpath = os.path.join(dat_dir, fname)
            
            if not os.path.exists(fpath):
                continue
            
            local_data = np.fromfile(fpath, dtype=np.float64)
            grid = local_data.reshape((lx, lz), order='F')
            
            global_data[ix*lx : (ix+1)*lx, iz*lz : (iz+1)*lz] = grid
            
    return global_data

def plot_physical_space(istep, var_name, npx, npz, lx, lz):
    # 1. 座標と時刻の読み込み
    x, z, t = load_coordinates()
    if x is None: return

    # 2. データの統合
    data = load_aggregated_step(istep, var_name, npx, npz, lx, lz)
    
    # 現在の時刻を取得 (istep番目の要素。0-based indexの場合は istep-1)
    current_time = t[istep] 

    # 3. 可視化
    plt.figure(figsize=(10, 7))
    
    # pcolormeshは座標(x, z)を直接扱えるため、非一様格子でも正しく表示可能
    # dataは (nx, nz) なので、プロット用に転置が必要な場合が多い
    # X, Z がセルセンターの場合、pcolormesh(X, Z, data.T) とします
    X, Z = np.meshgrid(x, z, indexing='ij')
    
    #pcm = plt.pcolormesh(X, Z, data, cmap='jet', shading='auto')
    pcm = plt.pcolormesh(X, Z, np.log10(data), cmap='jet', shading='auto')
    
    plt.rcParams['font.size']=11
    plt.rcParams['font.family']='Avenir'
    plt.rcParams['mathtext.fontset']='stix'
    plt.colorbar(pcm, label=f'Variable: {var_name}')
    plt.xlabel('Physical X [m]')
    plt.ylabel('Physical Z [m]')
    plt.ylim(0,1)
    plt.xlim(0,1)
    plt.title(f'Time = {current_time:.4f} s (Step: {istep:03d})')
    
    plt.tight_layout()
    plt.show()

# --- パラメータ設定 ---
if __name__ == "__main__":
    config = {
        "istep": 95,        # 表示したいステップ番号
        "var_name": "ro",   # 変数名
        "npx": 2, "npz": 2, # ランク分割数
        "lx": 200, "lz": 300 # 1ランクあたりの実格子数
    }
    
    plot_physical_space(**config)
