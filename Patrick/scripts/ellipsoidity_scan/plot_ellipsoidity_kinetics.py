import steps.interface
import steps.saving as stsave
from matplotlib import pyplot as plt
import seaborn as sns
import h5py
import numpy as np
import re
import math
import os
from os.path import join, dirname, abspath
import sys
# Add the project root directory to Python path
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../"))
sys.path.insert(0, project_root)

from Patrick.src.Utilities import get_repo_path

def get_ellipsoidity_from_mesh_name(mesh_path):
    """Extract ellipsoidity value from mesh filename."""
    filename = os.path.basename(mesh_path)
    match = re.search(r'ellipsoidity_(\d+\.?\d*)', filename)
    if match:
        return float(match.group(1))
    return None

def plot_ellipsoidity_kinetics(base_path, ellipsoidity_dir, figsize=(20, 15), exclude_species=None, output_filename='kinetics_all_species.png'):
    meshes_dir = join(base_path, "Patrick/meshes_ellipsoidity")

    print(f"Looking for results in: {ellipsoidity_dir}")
    print(f"Looking for meshes in: {meshes_dir}")
    
    # Get all mesh directories
    mesh_dirs = [d for d in os.listdir(ellipsoidity_dir) if d.startswith('mesh_')]
    mesh_dirs.sort(key=lambda x: int(x.split('_')[1]))  # Sort by mesh number
    print(f"Found {len(mesh_dirs)} mesh directories")
    
    # Store results for each metric
    ellipsoidity_values = []
    kinetic_data = {}
    
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
                    if species_name not in kinetic_data:
                        kinetic_data[species_name] = {}
                    
                    # Store all timepoints for all replicates
                    kinetic_data[species_name][ellipsoidity] = {
                        'time': res.time[0],
                        'data': res.data[:,:,0]  # Shape: (replicates, timepoints)
                    }
        except Exception as e:
            print(f"Error processing {result_path}: {str(e)}")
            continue
    
    if not kinetic_data:
        print("No kinetic data was collected. Check if the result files exist and contain the expected data.")
        return
        
    print(f"Collected data for {len(kinetic_data)} species")
    
    # Exclude specified species if requested
    if exclude_species is not None:
        kinetic_data = {k: v for k, v in kinetic_data.items() if k not in exclude_species}
        print(f"Plotting after excluding species: {exclude_species}. Remaining: {list(kinetic_data.keys())}")
    
    # Calculate grid layout
    num_species = len(kinetic_data)
    grid_size = math.ceil(math.sqrt(num_species))
    n_rows, n_cols = math.ceil(num_species / grid_size), grid_size
    
    # Create a single figure with subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()
    
    # Create a colormap for different ellipsoidity values
    # Use the first species to get ellipsoidity values (they should be the same for all species)
    first_species = list(kinetic_data.keys())[0]
    sorted_ellipsoidities = sorted(kinetic_data[first_species].keys())
    colors = plt.cm.viridis(np.linspace(0, 1, len(sorted_ellipsoidities)))
    
    # Plot each species in its own subplot
    for idx, (species_name, data) in enumerate(kinetic_data.items()):
        ax = axes[idx]
        
        for ellipsoidity, color in zip(sorted_ellipsoidities, colors):
            time = data[ellipsoidity]['time']
            all_replicates = data[ellipsoidity]['data']
            
            # Calculate mean and SEM across replicates
            mean_values = np.mean(all_replicates, axis=0)
            n_replicates = all_replicates.shape[0]
            sem_values = np.std(all_replicates, axis=0, ddof=1) / np.sqrt(n_replicates)
            
            # Plot mean line
            ax.plot(time, mean_values, 
                   label=f'Ellipsoidity = {ellipsoidity:.2f}',
                   color=color,
                   linewidth=2)
            
            # Plot standard error of the mean as shaded area
            ax.fill_between(time,
                          mean_values - sem_values,
                          mean_values + sem_values,
                          color=color,
                          alpha=0.2)
        
        ax.set_xlabel('Time [s]', fontsize=10)
        ax.set_ylabel(species_name, fontsize=10)
        # ax.set_title(species_name, fontsize=12)
        ax.grid(True, alpha=0.3)
    
    # Hide any unused subplots
    for i in range(idx + 1, len(axes)):
        fig.delaxes(axes[i])
    
    # Add a single legend for all lines, outside the grid to the right of the top row
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8, borderaxespad=0.)
    
    # Adjust layout to prevent overlap
    plt.tight_layout()
    
    # Save the plot
    output_path = join(ellipsoidity_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_path}")
    plt.close()

if __name__ == "__main__":
    base_path = get_repo_path()
    ellipsoidity_dir = join(base_path, "Patrick/saved_objects/ellipsoidity_few_slow_long")
    plot_ellipsoidity_kinetics(base_path, 
                               ellipsoidity_dir, 
                               figsize=(10, 5), 
                               exclude_species=["ERK", "EGF"], 
                               output_filename='kinetics_all_species.png') 