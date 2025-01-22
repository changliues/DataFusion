proximal_model_spec <- list(g.x = "g ~ x1+x2+x3+x1:x2+x1:x3+x2:x3+x1:x2:x3",
                            a.ex = "ae ~ x1e+x2e+x3e+x1e:x2e+x1e:x3e+x2e:x3e+x1e:x2e:x3e",
                            a.ox = "a ~ x1+x2+x3+x1:x2+x1:x3+x2:x3+x1:x2:x3",
                            eta = "w ~ x10e+x20e+x30e+x10e:x20e+x10e:x30e+x30e:x20e+x10e:x20e:x30e",
                            m.azo = "m ~ x1+x2+x3+z+a+I(z^2)+x1:x2+x1:x3+x1:z+x1:a+x1:I(z^2)+x2:x3+x2:z+x2:a+x2:I(z^2)+x3:z+x3:a+z:a+x3:I(z^2)+a:I(z^2)+x2:x3:z+x2:x3:a+x2:z:a+x2:x3:I(z^2)+x2:a:I(z^2)+x3:z:a+I(z^2):x3:a+
                            x1:x2:x3+x1:x3:z+x1:x3:a+x1:z:a+x1:x3:I(z^2)+x1:a:I(z^2)+x1:x2:z+x1:x2:a+x1:x2:I(z^2)+x1:x2:x3:z+x1:x2:x3:a+x1:x2:z:a+x1:x3:z:a+x2:x3:z:a+
                            x1:x2:x3:I(z^2)+x1:x2:a:I(z^2)+x1:I(z^2):x3:a+x2:I(z^2):x3:a+x1:x2:x3:z:a+x1:x2:I(z^2):x3:a",
                            y.azo = "y ~ x1+x2+x3+z+a+I(z^2)+x1:x2+x1:x3+x1:z+x1:a+x1:I(z^2)+x2:x3+x2:z+x2:a+x2:I(z^2)+x3:z+x3:a+z:a+x3:I(z^2)+a:I(z^2)+x2:x3:z+x2:x3:a+x2:z:a+x2:x3:I(z^2)+x2:a:I(z^2)+x3:z:a+I(z^2):x3:a+
                            x1:x2:x3+x1:x3:z+x1:x3:a+x1:z:a+x1:x3:I(z^2)+x1:a:I(z^2)+x1:x2:z+x1:x2:a+x1:x2:I(z^2)+x1:x2:x3:z+x1:x2:x3:a+x1:x2:z:a+x1:x3:z:a+x2:x3:z:a+
                            x1:x2:x3:I(z^2)+x1:x2:a:I(z^2)+x1:I(z^2):x3:a+x2:I(z^2):x3:a+x1:x2:x3:z:a+x1:x2:I(z^2):x3:a",
                            z.amo = "z ~ x1+x2+x3+m+a+I(m^2)+x1:x2+x1:x3+x1:m+x1:a+x1:I(m^2)+x2:x3+x2:m+x2:a+x2:I(m^2)+x3:m+x3:a+m:a+x3:I(m^2)+a:I(m^2)+x2:x3:m+x2:x3:a+x2:m:a+x2:x3:I(m^2)+x2:a:I(m^2)+x3:m:a+I(m^2):x3:a+
                            x1:x2:x3+x1:x3:m+x1:x3:a+x1:m:a+x1:x3:I(m^2)+x1:a:I(m^2)+x1:x2:m+x1:x2:a+x1:x2:I(m^2)+x1:x2:x3:m+x1:x2:x3:a+x1:x2:m:a+x1:x3:m:a+x2:x3:m:a+
                            x1:x2:x3:I(m^2)+x1:x2:a:I(m^2)+x1:I(m^2):x3:a+x2:I(m^2):x3:a+x1:x2:x3:m:a+x1:x2:I(m^2):x3:a",
                            m.ao = "m ~ x1+x2+x3+a+x1:x2+x1:x3+x1:a+x2:x3+x2:a+x3:a+x2:x3:a+x1:x2:x3+x1:x3:a+x1:x2:a+x1:x2:x3:a",
                            m.ae = "m ~ x1+x2+x3+a+x1:x2+x1:x3+x1:a+x2:x3+x2:a+x3:a+x2:x3:a+x1:x2:x3+x1:x3:a+x1:x2:a+x1:x2:x3:a")

saveRDS(proximal_model_spec, "proximal_model_spec_urban.rds")
