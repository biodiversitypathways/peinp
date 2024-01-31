wt_get_tasks <- function(project_id) {

  # User agent
  u <- getOption("HTTPUserAgent")
  if (is.null(u)) {
    u <- sprintf("R/%s; R (%s)",
                 getRversion(),
                 paste(getRversion(), R.version$platform, R.version$arch, R.version$os))
  }

  # Add wildRtrax version information:
  u <- paste0("wildRtrax ", as.character(packageVersion("wildRtrax")), "; ", u)

  # Modify the URL to include the dynamic project_id
  query_list = list(project_id = 2288)

  # Prepare temporary file:
  tmp <- tempfile(fileext = ".zip")
  # tmp directory
  td <- tempdir()

  # Create GET request
  taskz <- httr::POST(
    httr::modify_url("https://www-api.wildtrax.ca", path = "/bis/download-tasks-by-project"),
    query = query_list,
    accept = "application/json",
    httr::add_headers(Authorization = paste("Bearer", my_tok)),
    httr::user_agent(u),
    httr::write_disk(tmp),
    httr::progress()
  )

  if (taskz$status_code == 200) {
    # Success, continue processing
  } else {
    stop("The species table could not be downloaded.")
  }

  # Unzip
  unzip(tmp, exdir = td)

  # Remove abstract file
  tasks_table <- list.files(td, full.names = TRUE, recursive = TRUE)

  return(tasks_table)
}

wt_get_tasks(2288)
