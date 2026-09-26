# Small persistent config file for cross-session choices (e.g. a
# user-supplied cutadapt path), independent of any project directory.

app_config_path <- function() {
  file.path(path.expand("~"), ".arbusca", "config.yaml")
}

# Settings saved before the app was renamed (amfdada -> Arbusca) live in
# ~/.amfdada/; carry them over once so nobody has to re-enter them.
legacy_app_config_path <- function() {
  file.path(path.expand("~"), ".amfdada", "config.yaml")
}

read_app_config <- function() {
  path <- app_config_path()
  if (!file.exists(path) && file.exists(legacy_app_config_path())) {
    fs::dir_create(fs::path_dir(path))
    file.copy(legacy_app_config_path(), path)
  }
  if (!file.exists(path)) return(list())
  tryCatch(yaml::read_yaml(path), error = function(e) list())
}

write_app_config <- function(cfg) {
  path <- app_config_path()
  fs::dir_create(fs::path_dir(path))
  yaml::write_yaml(cfg, path)
}
