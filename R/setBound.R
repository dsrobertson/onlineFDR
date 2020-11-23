setBound <- function(alg, alpha, N) {
  bound <- switch(alg,
                  LOND = 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  LORD = 0.07720838 * log(pmax(seq_len(N + 1), 2))/(seq_len(N + 1) * 
                                                                      exp(sqrt(log(seq_len(N + 1))))),
                  LORDdep = 0.139307 * alpha/(b0 * seq_len(N) * (log(pmax(seq_len(N), 2)))^3),
                  SAFFRON = 0.4374901658/(seq_len(N)^(1.6)),
                  ADDIS = 0.4374901658/(seq_len(N + 1)^(1.6)),
                  LONDSTAR = 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  LORDSTAR = 0.07720838 * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  SAFFRONSTAR = 0.4374901658/(seq_len(N + 1)^(1.6))
  )
}
