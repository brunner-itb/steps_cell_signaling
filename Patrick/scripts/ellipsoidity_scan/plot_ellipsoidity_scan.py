import steps.interface
import steps.saving as stsave
import os
import sys

# Add the project root directory to Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
sys.path.insert(0, project_root)

from Patrick.src.Utilities import get_repo_path
from matplotlib import pyplot as plt
import seaborn as sns
import h5py
import numpy as np
import re
import math
from os.path import join




def get_ellipsoidity_from_mesh_name(mesh_path):
    """Extract ellipsoidity value from mesh filename."""
    filename = os.path.basename(mesh_path)
    match = re.search(r'ellipsoidity_(\d+\.?\d*)', filename)
    if match:
        return float(match.group(1))
    return None

def plot_ellipsoidity_metrics(base_path, ellipsoidity_dir):
    meshes_dir = join(base_path, "Patrick/meshes_ellipsoidity")

    print(f"Looking for results in: {ellipsoidity_dir}")
    print(f"Looking for meshes in: {meshes_dir}")
    
    # Get all mesh directories
    mesh_dirs = [d for d in os.listdir(ellipsoidity_dir) if d.startswith('mesh_')]
    mesh_dirs.sort(key=lambda x: int(x.split('_')[1]))  # Sort by mesh number
    print(f"Found {len(mesh_dirs)} mesh directories")
    
    # Store results for each metric
    ellipsoidity_values = []
    metric_data = {}
    
    # Process each mesh result
    for mesh_dir in mesh_dirs:
        result_path = join(ellipsoidity_dir, mesh_dir, 'result')
        if not os.path.exists(result_path + '.h5'):
            print(f"No result file found for {mesh_dir}")
            continue
            
        # Get ellipsoidity value from mesh name
        mesh_idx = mesh_dir.split('_')[1]
        mesh_files = [f for f in os.listdir(meshes_dir) if f.startswith('ellipsoidity_') and f.endswith('.inp')]
        mesh_files.sort(key=lambda x: float(re.search(r'ellipsoidity_(\d+\.?\d*)', x).group(1)))
        
        if int(mesh_idx) >= len(mesh_files):
            print(f"No matching mesh file found for {mesh_dir}")
            continue
            
        mesh_path = join(meshes_dir, mesh_files[int(mesh_idx)])
        ellipsoidity = get_ellipsoidity_from_mesh_name(mesh_path)
        if ellipsoidity is None:
            print(f"Could not extract ellipsoidity from {mesh_path}")
            continue
            
        print(f"Processing mesh {mesh_dir} with ellipsoidity {ellipsoidity}")
        ellipsoidity_values.append(ellipsoidity)
        
        # Load results
        try:
            with stsave.HDF5Handler(result_path) as hdf:
                results = hdf["ellipsoidity"].results
                
                # Process each species/metric
                for res in results:
                    species_name = re.search(r'\.(.*?)\.', res.labels[0]).group(1)
                    if species_name not in metric_data:
                        metric_data[species_name] = []
                    
                    # Store all replicate values
                    replicate_values = res.data[:,:,0]  # Shape: (replicates, timepoints)
                    final_values = replicate_values[:, -1]  # Get last timepoint for each replicate
                    metric_data[species_name].append(final_values)
        except Exception as e:
            print(f"Error processing {result_path}: {str(e)}")
            continue
    
    if not metric_data:
        print("No metric data was collected. Check if the result files exist and contain the expected data.")
        return
        
    print(f"Collected data for {len(metric_data)} metrics")
    
    # Sort data by ellipsoidity
    sort_idx = np.argsort(ellipsoidity_values)
    ellipsoidity_values = np.array(ellipsoidity_values)[sort_idx]
    for species in metric_data:
        metric_data[species] = np.array(metric_data[species])[sort_idx]
    
    # Create plots
    num_metrics = len(metric_data)
    grid_size = math.ceil(math.sqrt(num_metrics))
    n_rows, n_cols = grid_size, math.ceil(num_metrics / grid_size)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, (species, values) in enumerate(metric_data.items()):
        ax = axes[idx]
        
        # Calculate mean and std across replicates
        mean_values = np.mean(values, axis=1)
        std_values = np.std(values, axis=1)
        
        # Plot mean with error bars
        ax.errorbar(ellipsoidity_values, mean_values, yerr=std_values, 
                   fmt='o-', capsize=5, capthick=1, elinewidth=1,
                   label=species)
        
        ax.set_xlabel('Ellipsoidity')
        ax.set_ylabel('Final Concentration')
        ax.set_title(species)
        ax.grid(True)
    
    # Hide unused subplots
    for i in range(idx + 1, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    output_path = join(ellipsoidity_dir, 'ellipsoidity_metrics.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    plt.close()

if __name__ == "__main__":
    base_path = get_repo_path()
    ellipsoidity_dir = join(base_path, "Patrick/saved_objects/ellipsoidity")
    plot_ellipsoidity_metrics(base_path, ellipsoidity_dir) 