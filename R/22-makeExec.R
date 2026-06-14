# Make python executable for standard (feedforward) implementation

#' @title Generate Python Executable for Neural Network
#' @description
#' Creates a Python executable for feedforward neural network implementation
#' by combining multiple script parts into a complete TensorFlow model.
#'
#' @param tmpDir Temporary directory for script generation
#' @param FFdepth Number of hidden layers in the feedforward network
#' @param model_type Type of neural network model. Must be one of:
#' 'BlankClass', 'ClassBlank', 'BlankCox', 'ClassClass', or 'ClassCox'
#'
#' @return None. Generates Python scripts in the specified temporary directory.
#'
#' @examples
#' \dontrun{
#' # Create a 2-layer classifier model
#' makeExec("/tmp/model/", 2, "ClassClass")
#' }
#'
#' @export
makeExec <- function(
  # tmpDir,
  FFdepth,
  model_type = c(
    'ClassClass',
    'ClassCox',
    'ClassBlank',
    'BlankClass',
    'BlankCox'
  ),
  DEGAS.toolsPath = system.file(
    "DEGAS_tools",
    package = "DEGAS",
    mustWork = TRUE
  ),
  architecture = c(
    "Standard", # ! makeExec
    "DenseNet" # ! makeExec2
  ),
  ...
) {
  model_type <- match.arg(model_type)
  architecture <- match.arg(architecture)

  p1_p3_filenames <- paste0(model_type, c('MTL_p1.py', 'MTL_p3.py'))

  #   system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p1.py ', tmpDir))
  p1_content <- read_py_as_chr(p1_p3_filenames[1])
  #   system(paste0('cp ', DEGAS.toolsPath, model_type, 'MTL_p3.py ', tmpDir))
  p3_content <- read_py_as_chr(p1_p3_filenames[2])

  p2_content <- generateMTLp2(FFdepth = FFdepth, architecture = architecture)

  #   fout = file(paste0(tmpDir, model_type, 'MTL_p2.py'))
  #   writeLines(outlines, fout)
  #   close(fout)

  p4_content <- generateMTLp4(FFdepth = FFdepth, model_type = model_type)

  #   fout = file(paste0(tmpDir, model_type, 'MTL_p4.py'))
  #   writeLines(outlines, fout)
  #   close(fout)

  out_filename <- paste0(model_type, 'MTL.py')
  out_content <- paste(
    p1_content,
    p2_content,
    p3_content,
    p4_content,
    sep = '\n',
    collapse = '\n'
  )

  # Prepend graph reset and config for safe inline execution via reticulate
  out_content <- paste(
    '# Reset TF graph to avoid variable accumulation in reticulate session',
    'tf.reset_default_graph()',
    'config = tf.ConfigProto()',
    out_content,
    sep = '\n',
    collapse = '\n'
  )

  #   writeLines(out_content, con = file.path(DEGAS.toolsPath, out_filename))

  out_content
}
