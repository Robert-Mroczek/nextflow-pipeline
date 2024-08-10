process dada2_varcall { 
    cpus 3
    publishDir("./results/dadavariantcall", mode: 'copy')
    input:
        val(fow_only)
    output: 
        tuple path('*_sequence_table'), path('*.pdf')
    script: 
        sample_map = [:]
        sample_csv = file("./sample_setup.csv") 
        allLines = sample_csv.readLines()
        for (line : allLines){
            splited_line = line.split(",") 
            sample_map[splited_line[0]] = splited_line[1] 
        }
        dada_sample_list = []
        fow_only.eachWithIndex { item, index ->
            if (item[1]){
                dada_sample_list.add(item[0])
                }
                
        }
        dada_sample_list.eachWithIndex{ item, index ->
            if (sample_map[item[0]] == 'sample1'){dada_sample1 = item[1]} 
            if (sample_map[item[0]] == 'sample2'){dada_sample2 = item[1]} 
            if (sample_map[item[0]] == 'control'){dada_control = item[1]} 
            }

        """
        Rscript "$baseDir/templates/dada.R" $dada_sample1 $dada_sample2 $dada_control 
        """
}