Athey_model_spec <- list(y.m0o = "y ~ x1*x2*x3*x4*m*I(x3^2)*I(m^2)")
naive_model_spec <- list(y.0o = "y ~ x1*x2*x3*x4*I(x3^2)",
                         a.e = "a ~ x1*x2*x3*x4*I(x3^2)",
                         a.o = "a ~ x1*x2*x3*x4*I(x3^2)",
                         g = "gn ~ x1*x2*x3*x4*I(x3^2)",
                         m.0e = "m ~ x1*x2*x3*x4*I(x3^2)",
                         m.0o = "m ~ x1*x2*x3*x4*I(x3^2)")


saveRDS(Athey_model_spec, "Athey_model_spec.rds")
saveRDS(naive_model_spec, "naive_model_spec.rds")
