if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")

setRepositories()

install.packages("rsconnect")


"/home/kvalem/R/x86_64-pc-linux-gnu-library/4.4/pavian/shinyapp"


rsconnect::setAccountInfo(name='erikakvalem',
                          token='36DDE2F4FDD1FF192823A1A7208DA514',
                          secret='3Ll/Z494CxtDJ3IJUctxFm3TZVxo6HIbbiu7D2wx')
rsconnect::deployApp("/home/kvalem/R/x86_64-pc-linux-gnu-library/4.4/pavian/shinyapp")
