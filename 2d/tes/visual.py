import numpy as np
import matplotlib.pyplot as plt
import os

def aggregate_and_plot_2d(base_fname, npx, npz, lx, lz):
    """
    全ランクのバイナリファイルを読み込み、1つのグローバル配列に統合して可視化する
    
    Parameters:
    -----------
    base_fname : str
        ファイル名のプレフィックス (例: 'data')
    npx, npz   : int
        X方向およびZ方向の並列分割数
    lx, lz     : int
        各ランクが持つ実計算領域のサイズ (lix-2*mg, liz-2*mg)
    """
    
    # 1. グローバル配列の確保 (全領域をカバーするサイズ)
    global_data = np.zeros((npx * lx, npz * lz), dtype=np.float64)
    
    print(f"Global array initialized: {global_data.shape}")

    # 2. 各ランクのファイルをループして読み込み
    for i in range(npx):
        for j in range(npz):
            # Fortranの write(rank_fname,'(A,"_",i0,"_",i0,".dat")') に合わせる
            file_path = f"{base_fname}_{i}_{j}.dat"
            
            if not os.path.exists(file_path):
                print(f"Warning: {file_path} not found. Skipping...")
                continue
            
            # バイナリ読み込み
            local_data = np.fromfile(file_path, dtype=np.float64)
            
            # 形状を復元 (Fortran order)
            local_grid = local_data.reshape((lx, lz), order='F')
            
            # グローバル配列の該当箇所に挿入
            x_start, x_end = i * lx, (i + 1) * lx
            z_start, z_end = j * lz, (j + 1) * lz
            global_data[x_start:x_end, z_start:z_end] = local_grid

    # 3. 可視化
    plt.figure(figsize=(10, 8))
    plt.rcParams['font.size']=11
    plt.rcParams['font.family']='Avenir'
    plt.rcParams['mathtext.fontset']='stixsans'
    
    # Fortranのインデックス(x, z)を画像(row, col)に変換するために転置 (.T)
    # origin='lower' で (0,0) を左下に配置
    img = plt.imshow(global_data.T, origin='lower', cmap='jet', aspect='equal',
                     extent=[0, npx*lx, 0, npz*lz])
    
    plt.colorbar(img, label='Physical Value')
    plt.xlabel('Global X Index')
    plt.ylabel('Global Z Index')
    plt.title(f'Aggregated Simulation Result ({npx}x{npz} ranks)')
    
    plt.tight_layout()
    plt.show()

# --- 設定値 (適宜書き換えてください) ---
if __name__ == "__main__":
    params = {
        "base_fname": "dat/ro",  # Fortranのfnameに相当
        "npx": 3,              # X方向のランク数
        "npz": 2,              # Z方向のランク数
        "lx": 48,             # 1ランクあたりのXサイズ (lix - 2*mg)
        "lz": 50              # 1ランクあたりのZサイズ (liz - 2*mg)
    }
    
    aggregate_and_plot_2d(**params)
