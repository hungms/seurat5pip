#' integrate_by_seurat
#' 
#' @param obj Seurat object
#' @param split.by Split by
#' @param vars.to.regress Variables to regress out
#' @param method Method for integration
#' @param normalization.method Normalization method
#' @param return.only.var.genes Return only variable genes
#' @param k.weight K weight
#' 
#' @return Seurat object
int_by_seurat <- function(obj, split.by = NULL, normalization.method = "SCT", vars.to.regress = NULL, method = "rpca", return.only.var.genes = T, k.weight = 100){

    split.by <- validate_split.by(split.by, obj)

    options(Seurat.object.assay.version = "v5")
    features <- VariableFeatures(obj)
    obj.list <- SplitObject(obj, split.by = split.by)

    if(normalization.method == "SCT"){
        obj.list <- lapply(X = obj.list, function(x) SCTransform(x, vars.to.regress = vars.to.regress, verbose = FALSE))
        obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)}

    obj.list <- lapply(X = obj.list, FUN = RunPCA, features = features)
    anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = normalization.method, anchor.features = features, reduction = method)
    if(return.only.var.genes){
        integrated <- IntegrateData(anchorset = anchors, normalization.method = normalization.method, features.to.integrate = features, k.weight = k.weight)}
    else{
        integrated <- IntegrateData(anchorset = anchors, normalization.method = normalization.method, features.to.integrate = rownames(obj), k.weight = k.weight)}
    VariableFeatures(integrated[["SCT"]]) <- features
    #options(Seurat.object.assay.version = "v5")
    return(integrated)}
