import sys
import argparse
import cnmf
import yaml
import os
import pandas as pd



# Change path to wherever you have repo locally
sys.path.append('/project/GCRB/Hon_lab/s223695/Data_project/jamboree_2025/cNMF_benchmarking/cNMF_benchmarking_pipeline')


from Inference.src import (
    run_cnmf_consensus, get_top_indices_fast, annotate_genes_to_excel, \
    rename_and_move_files_NMF, rename_all_NMF, compile_results
)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--counts_fn', type=str, required=True)  
    parser.add_argument('--output_directory', type=str, required=True)
    parser.add_argument('--run_name', type=str, required=True)

    parser.add_argument('--K', nargs='*', type=int, default=None)
    parser.add_argument('--numiter', type = int, default = 10)
    parser.add_argument('--densify', action='store_true') # Flag: include --densify for True, omit for False
    parser.add_argument('--tpm_fn', type = str, default = None)
    parser.add_argument('--numhvgenes', type = int, default = 5451)
    parser.add_argument('--genes_file', type = str, default = None)
    parser.add_argument('--init', type = str, default = 'random')
    parser.add_argument('--loss', default = 'frobenius')
    parser.add_argument('--algo', type = str, default = 'mu')
    parser.add_argument('--mode', type = str, default = 'batch')
    parser.add_argument('--tol', type = float, default = 1e-4)
    parser.add_argument('--n_jobs', type = int, default = -1)
    parser.add_argument('--seed', type = int, default = 14)
    parser.add_argument('--use_gpu', action="store_true") # Flag: include --use_gpu for True, omit for False
    parser.add_argument('--alpha_usage', type = float, default = 0.0)
    parser.add_argument('--alpha_spectra', type = float, default = 0.0)
    parser.add_argument('--l1_ratio_usage', type = float, default = 0.0)
    parser.add_argument('--l1_ratio_spectra', type = float, default = 0.0)
    parser.add_argument('--online_usage_tol', type = float, default = 0.05)
    parser.add_argument('--online_spectra_tol', type = float, default = 0.05)
    parser.add_argument('--fp_precision', type = str, default = 'float')
    parser.add_argument('--batch_max_iter', type = int, default = 500)
    parser.add_argument('--batch_hals_tol', type = float, default = 0.05)
    parser.add_argument('--batch_hals_max_iter', type = int, default = 200)
    parser.add_argument('--online_max_pass', type = int, default = 20)
    parser.add_argument('--online_chunk_size', type = int, default = 5000)
    parser.add_argument('--online_chunk_max_iter', type = int, default = 200)
    parser.add_argument('--shuffle_cells', action="store_true")
    parser.add_argument('--sk_cd_refit', action="store_true")
    
    parser.add_argument('--parallel_running', action="store_true")
    parser.add_argument('--num_gene', type = int, default = 300)

    parser.add_argument('--sel_thresh', nargs='*', type=float, default=None)  

    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        k_value = [30, 50, 60, 80, 100, 200, 250, 300]
    else:
        k_value = args.K

    if args.sel_thresh is None:
        sel_thresh_value = [0.2, 2.0]
    else:
        sel_thresh_value = args.sel_thresh


    # save comfigs used         
    os.makedirs((f'{args.output_directory}/{args.run_name}'), exist_ok=True)
    args_dict = vars(args)
    with open(f'{args.output_directory}/{args.run_name}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)

    # running 
    cnmf_obj = cnmf.cNMF(output_dir=args.output_directory, name=args.run_name)
    print(f"{args.counts_fn} used for inference" )


    cnmf_obj.prepare(counts_fn=args.counts_fn, components=k_value, n_iter=args.numiter, densify=args.densify, tpm_fn=args.tpm_fn, num_highvar_genes=args.numhvgenes, genes_file=args.genes_file,
                init = args.init,  beta_loss = args.loss, 
                algo = args.algo, mode = args.mode, tol=args.tol, n_jobs=args.n_jobs, 
                seed = args.seed,  use_gpu = args.use_gpu, 
                alpha_usage = args.alpha_usage, alpha_spectra = args.alpha_spectra, 
                l1_ratio_usage = args.l1_ratio_usage, l1_ratio_spectra = args.l1_ratio_spectra,
                online_usage_tol = args.online_usage_tol, online_spectra_tol = args.online_spectra_tol,
                fp_precision = args.fp_precision, 
                batch_max_iter = args.batch_max_iter, batch_hals_tol = args.batch_hals_tol, batch_hals_max_iter = args.batch_hals_max_iter,
                online_max_pass = args.online_max_pass, online_chunk_size = args.online_chunk_size, online_chunk_max_iter = args.online_chunk_max_iter,
                shuffle_cells = args.shuffle_cells, sk_cd_refit =args.sk_cd_refit )



    cnmf_obj.factorize(total_workers=1)

    cnmf_obj.combine()

    cnmf_obj.k_selection_plot()

    # Consensus plots with all k to choose thresh
    run_cnmf_consensus(cnmf_obj, 
                       components=k_value, 
                       density_thresholds=sel_thresh_value)

#     # Save all cNMF scores in separate mudata objects
#     compile_results(args.output_directory, args.run_name, components= k_value, sel_thresh = sel_thresh_value )

#     # annotation for all K
#     os.makedirs((f'{args.output_directory}/{args.run_name}/Annotation'), exist_ok=True)

#     for i in sel_thresh_value:
#         for k in args.K:
#             df = pd.read_csv('{output_directory}/{run_name}/{run_name}.gene_spectra_scores.k_{k}.dt_{sel_thresh}.txt'.format(
#                                                                                     output_directory=args.output_directory,
#                                                                                     run_name = args.run_name,
#                                                                                     k=k,
#                                                                                     sel_thresh = str(i).replace('.','_')),
#                                                                                     sep='\t', index_col=0)   
#             overlap = get_top_indices_fast(df, gene_num=args.num_gene)
#             annotate_genes_to_excel(overlap, f'{args.output_directory}/{args.run_name}/Annotation/{k}_{i}.xlsx')



#     # combine the parallel run K value into "run_name_all" file 
#     if args.parallel_running and isinstance(k_value, int):

#         file_name_input_new = f"{file_name_input}_{k_value}.spectra.k_{k_value}.iter"
#         file_name_output_new = f"{args.run_name}_all.spectra.k_{k_value}.iter"

#         source_folder_new = f"{source_folder}_{k_value}/cnmf_tmp"

#         rename_all_NMF(source_folder = args.output_directory ,
#                                 destination_folder = f"{args.output_directory}_all/cnmf_tmp",
#                                 file_name_input = args.run_name,
#                                 file_name_output = f"{args.run_name}_all",
#                                 len = args.numiter, 
#                                 components = k_value)



