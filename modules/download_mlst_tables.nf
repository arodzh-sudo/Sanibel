process download_mlst_tables {
    storeDir "${projectDir}/db/mlst_tables"

    output:
        path 'neisseria.txt',   emit: neisseria
        path 'hinfluenzae.txt', emit: hinfluenzae

    script:
    """
    wget -q -O neisseria.txt   https://raw.githubusercontent.com/tseemann/mlst/master/db/pubmlst/neisseria/neisseria.txt
    wget -q -O hinfluenzae.txt https://raw.githubusercontent.com/tseemann/mlst/master/db/pubmlst/hinfluenzae/hinfluenzae.txt
    """
}
