# VLSM: use chi squared between AD and CN patients
using StatsModels, CSV, DataFrames, GLM, NIfTI, ArgParse
using MultipleTesting, HypothesisTests


function run_test(df, xvar, outcome, covars, method)
    # `df`: dataframe, must contain outcome and the covariates
    # `xvar`: a vector, same length as the number of rows in `df`
    df.x = xvar
    no_na = dropmissing(df, [outcome])
    no_na.y = no_na[:, outcome]

    all_vars = [Term(Symbol("x"))]
    for i in 1:length(covars)
        push!(all_vars, Term(Symbol(covars[i])))
    end

    if (method == "chisq")
        # create a 2x2 contingency table
        contig = [0 0; 0 0]
        contig[1, 1] = sum((no_na.x .== 0).*(no_na.y .== 0))
        contig[1, 2] = sum((no_na.x .== 0).*(no_na.y .== 1))
        contig[2, 1] = sum((no_na.x .== 1).*(no_na.y .== 0))
        contig[2, 2] = sum((no_na.x .== 1).*(no_na.y .== 1))
        model = ChisqTest(contig)
        pval = pvalue(model)
        statval = model.stat
    else
        model = lm(FormulaTerm(Term(Symbol("y")), Tuple(all_vars)), no_na)
        ctab = coeftable(model)
        ctab = DataFrame(ctab)
        pval = ctab[(ctab.Name .== "x"), "Pr(>|t|)"][1]
        statval = ctab[(ctab.Name .== "x"), "Coef."][1]
    end

    # create a vector of voxels and covariates
    var_stat = zeros(2)
    var_stat[1] = pval
    var_stat[2] = statval
    return var_stat
end


function sum_map(all_files, threshold, save_path) 
    # find out which voxels happen at a frequency above the threshold
    one_img = niread(all_files[1])
    sum_img = one_img.raw
    for i in 2:length(all_files)
        img = niread(all_files[i])
        sum_img = sum_img + img.raw
    end
    save_nii(sum_img, one_img, save_path)
    proceed_map = sum_img .> threshold*length(all_files)
    return proceed_map
end


function obtain_coord(proceed_map)
    dimen = size(proceed_map)
    N = floor(Int, sum(proceed_map))
    coord = zeros(Int, (N, 3))
    index = 1
    for k in 1:dimen[3]
        for j in 1:dimen[2]
            for i in 1:dimen[1]
                if (proceed_map[i, j, k] == true)
                    coord[index, 1] = i
                    coord[index, 2] = j
                    coord[index, 3] = k
                    index = index + 1
                end
            end
        end
    end
    return coord
end

function obtain_matrix(all_files, coord)
    N = length(all_files)
    L = size(coord)[1]
    all_mat = zeros(N, L)
    for i in 1:N
        img = niread(all_files[i])
        for j in 1:L
            x = coord[j, 1]
            y = coord[j, 2]
            z = coord[j, 3]
            all_mat[i, j] = img.raw[x, y, z]
        end
    end
    return all_mat
end

function run_test_all(all_df, all_mat, outcome, covar, method)
    N = size(all_mat)[1]
    L = size(all_mat)[2]
    M = length(covar)
    out_stat = zeros(L, 2)
    for i in 1:L
        out = run_test(all_df, all_mat[:, i], outcome, covar, method)
        out_stat[i, :] = out
    end

    out_stat[:, 1] = adjust(out_stat[:, 1], BenjaminiHochberg())
    return out_stat
end


function coord_to_map(vec, coord, img_size)
    # convert tables of p values or other statistical outcomes into maps
    N = size(vec)[1]
    TY = size(vec)[2]
    img = ones(img_size[1], img_size[2], img_size[3], TY)
    for j in 1:TY
        for i in 1:N
            x = coord[i, 1]
            y = coord[i, 2]
            z = coord[i, 3]
            img[x, y, z, j] = vec[i, j]
        end
    end
    return img
end

function save_nii(arr, nifti_ob, save_path)
    arr_nib = NIVolume(nifti_ob.header, arr)
    niwrite(save_path, arr_nib)
end

function vox_wise_test(df_path, outcome, save_dir, covar, method, vox_thres)
    # save_dir: save 2 images: stat_img.nii.gz (map of hazard ratio),
    # p_val.nii.gz (p value map after FDR correction)
    println("reading csv")
    all_df = CSV.read(df_path, DataFrame)
    all_paths = all_df.paths

    println("loading maps")
    proceed_map = sum_map(all_paths, vox_thres, save_dir*"/freq_map.nii.gz")
    coord_one = obtain_coord(proceed_map);
    all_mat = obtain_matrix(all_paths, coord_one);

    println("voxel wise statistical testing")
    out_mat = run_test_all(all_df, all_mat, outcome, covar, method)

    println("saving output")
    one_img = niread(all_paths[1])
    out_maps = coord_to_map(out_mat, coord_one, size(one_img))

    M = length(covar)
    save_nii(out_maps[:, :, :, 1], one_img, save_dir*"/pval_img.nii.gz")
    save_nii(out_maps[:, :, :, 2], one_img, save_dir*"/stat_img.nii.gz")
    println("timing analysis")
end


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--dfpath"
            help = "path to the image and behavior dataframes"
            required = true
        "--outcome"
            help = "which outcome to predict"
            required = true
        "--save_dir"
            help = "where to save the data"
            required = true
        "--covar"
            help = "covariates (comma separated)"
            default = "Age,Sex"
        "--test"
            help = "statistical test"
            default = "lingress"
    end
    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    if !isdir(parsed_args["save_dir"])
        mkdir(parsed_args["save_dir"])
    end
    covar = collect(eachsplit(parsed_args["covar"], ","))
    vox_wise_test(parsed_args["dfpath"], parsed_args["outcome"],
                  parsed_args["save_dir"], covar, parsed_args["test"], 0.025)
end

@time main()
