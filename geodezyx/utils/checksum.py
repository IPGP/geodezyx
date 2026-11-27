import os
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm


def calculate_sha256(filepath):
    """Calcule le hash SHA-256 d'un fichier."""
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(8192):
            sha256.update(chunk)
    return sha256.hexdigest()


def _scandir_recursive(root_dir):
    """
    Générateur efficace utilisant os.scandir() - cache les métadonnées.
    Meilleur pour millions de fichiers: ~2-3x plus rapide que os.walk.
    """
    try:
        for entry in os.scandir(root_dir):
            if entry.is_file(follow_symlinks=False):
                yield entry.path
            elif entry.is_dir(follow_symlinks=False):
                yield from _scandir_recursive(entry.path)
    except (PermissionError, OSError):
        pass


def _count_files_scandir(root_dir):
    """Compte les fichiers avec scandir (cache les métadonnées)."""
    count = 0
    try:
        for entry in os.scandir(root_dir):
            if entry.is_file(follow_symlinks=False):
                count += 1
            elif entry.is_dir(follow_symlinks=False):
                count += _count_files_scandir(entry.path)
    except (PermissionError, OSError):
        pass
    return count


def get_files_hashes(root_dir, output_summary_file=None, flush_interval=1000, use_scandir=True):
    """
    Calcule les hashes SHA-256 de tous les fichiers.
    
    Args:
        root_dir: Répertoire racine
        output_summary_file: Fichier de sortie (optionnel)
        flush_interval: Écrire tous les N fichiers
        use_scandir: True → os.scandir (recommandé pour millions), False → os.walk
    """
    from concurrent.futures import as_completed
    
    outdic = {}
    count = 0
    
    # Choisir la stratégie de scan
    if use_scandir:
        total_files = _count_files_scandir(root_dir)
        file_generator = _scandir_recursive(root_dir)
    else:
        walk_data = list(os.walk(root_dir))
        total_files = sum(len(files) for _, _, files in walk_data)
        file_generator = (
            os.path.join(dirpath, filename)
            for dirpath, _, filenames in walk_data
            for filename in filenames
        )
    
    with ThreadPoolExecutor() as executor:
        # Soumettre toutes les tâches
        futures_to_path = {}
        for filepath in file_generator:
            future = executor.submit(calculate_sha256, filepath)
            futures_to_path[future] = filepath
        
        # Traiter les résultats au fur et à mesure
        with tqdm(total=total_files, desc="Hashing files", unit="file") as pbar:
            for future in as_completed(futures_to_path):
                filepath = futures_to_path[future]
                abs_path = os.path.abspath(filepath)
                file_hash = future.result()
                rel_path = os.path.relpath(filepath, root_dir)
                outdic[rel_path] = file_hash
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
