.onLoad <- function(libname, pkgname) {
  # 1) absoluten Pfad auf die nlopt.dll ermitteln
  dll <- system.file("bin", .Platform$r_arch, "nlopt.dll", package = pkgname)
  if (file.exists(dll)) {
    # 2) VOR dem Laden von QuoteDynamics.dll einmal dyn.load aufrufen
    dyn.load(dll, local = FALSE, now = TRUE)
  }
}
