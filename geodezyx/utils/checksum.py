import os
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm


def calculate_sha256(filepath):
    """Calcule le hash SHA-256 d'un fichier."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):  # Lecture par blocs pour optimiser la mémoire
            sha256.update(chunk)
    return sha256.hexdigest()

def get_files_hashes(root_dir, output_summary_file=None, flush_interval=1000):
    from concurrent.futures import as_completed
    
    hashes = {}
    outdic = {}
    count = 0
    
    # Stocker os.walk une seule fois
    walk_data = list(os.walk(root_dir))
    total_files = sum(len(files) for _, _, files in walk_data)
    
    with ThreadPoolExecutor() as executor:
        # Soumettre toutes les tâches
        for dirpath, _, filenames in walk_data:
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                rel_path = os.path.relpath(filepath, root_dir)
                hashes[rel_path] = (executor.submit(calculate_sha256, filepath), filepath)
        
        # Traiter les résultats au fur et à mesure
        with tqdm(total=total_files, desc="Hashing files", unit="file") as pbar:
            futures_to_path = {hashes[k][0]: k for k in hashes}
            
            for future in as_completed(futures_to_path):
                rel_path = futures_to_path[future]
                abs_path = hashes[rel_path][1]
                outdic[rel_path] = future.result()
                pbar.update(1)
                count += 1
                
                # Écrire périodiquement
                if output_summary_file and count % flush_interval == 0:
                    _write_batch(outdic, output_summary_file, count == flush_interval)
    
    # Écrire les résultats restants
    if output_summary_file and outdic:
        _write_batch(outdic, output_summary_file, False)
    
    return outdic


def _write_batch(results, output_file, first_write=False):
    """Écrit un batch de résultats dans le fichier."""
    if not results:
        return
    
    mode = "w" if first_write else "a"
    with open(output_file, mode) as f:
        for rel_path, file_hash in results.items():
            abs_path = os.path.abspath(rel_path)
            f.write(f"{file_hash} {abs_path}\n")
    
    results.clear()


# def compare_directories(dir1, dir2):
#     """Compare deux arborescences et retourne les différences."""
#     hashes1 = get_files_hashes(dir1)
#     hashes2 = get_files_hashes(dir2)
#
#     # Clés uniques à chaque arborescence
#     only_in_dir1 = set(hashes1.keys()) - set(hashes2.keys())
#     only_in_dir2 = set(hashes2.keys()) - set(hashes1.keys())
#
#     # Fichiers communs mais avec des hashs différents
#     common_files = set(hashes1.keys()) & set(hashes2.keys())
#     different_hashes = {
#         file: (hashes1[file], hashes2[file])
#         for file in common_files
#         if hashes1[file] != hashes2[file]
#     }
#
#     return {
#         "fichiers_uniquement_dans_dir1": list(only_in_dir1),
#         "fichiers_uniquement_dans_dir2": list(only_in_dir2),
#         "fichiers_differents": different_hashes,
#     }
#
# # Exemple d'utilisation
# if __name__ == "__main__":
#     dir1 = "/chemin/vers/arborescence1"
#     dir2 = "/chemin/vers/arborescence2"
#
#     differences = compare_directories(dir1, dir2)
#
#     print("=== Fichiers uniquement dans dir1 ===")
#     for file in differences["fichiers_uniquement_dans_dir1"]:
#         print(f"- {file}")
#
#     print("\n=== Fichiers uniquement dans dir2 ===")
#     for file in differences["fichiers_uniquement_dans_dir2"]:
#         print(f"- {file}")
#
#     print("\n=== Fichiers avec des contenus différents ===")
#     for file, (hash1, hash2) in differences["fichiers_differents"].items():
#         print(f"- {file} (hash_dir1: {hash1}, hash_dir2: {hash2})")
