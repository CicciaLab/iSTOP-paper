# ---- Figure 4 ----
# The remaining species datasets are necessary to reproduce Figure 4
# Worm ~ 450 MB
read_csv('data/CDS/Celegans-ce11.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Celegans.UCSC.ce11::Celegans, cores = cores) %>%
  locate_PAM(BSgenome.Celegans.UCSC.ce11::Celegans) %>%
  write_csv('data/iSTOP/Celegans-ce11.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Celegans-ce11.csv')

# Fly ~ 650 MB
read_csv('data/CDS/Dmelanogaster-dm6.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Dmelanogaster.UCSC.dm6::Dmelanogaster, cores = cores) %>%
  locate_PAM(BSgenome.Dmelanogaster.UCSC.dm6::Dmelanogaster) %>%
  write_csv('data/iSTOP/Dmelanogaster-dm6.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Dmelanogaster-dm6.csv')

# Fish ~ 580 MB
read_csv('data/CDS/Drerio-danRer10.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Drerio.UCSC.danRer10::Drerio, cores = cores) %>%
  locate_PAM(BSgenome.Drerio.UCSC.danRer10::Drerio) %>%
  write_csv('data/iSTOP/Drerio-danRer10.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Drerio-danRer10.csv')

# Mouse ~ 720 MB
read_csv('data/CDS/Mmusculus-mm10.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, cores = cores) %>%
  locate_PAM(BSgenome.Mmusculus.UCSC.mm10::Mmusculus) %>%
  write_csv('data/iSTOP/Mmusculus-mm10.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Mmusculus-mm10.csv')

# Rat ~ 420 MB
read_csv('data/CDS/Rnorvegicus-rn6.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Rnorvegicus.UCSC.rn6::Rnorvegicus, cores = cores) %>%
  locate_PAM(BSgenome.Rnorvegicus.UCSC.rn6::Rnorvegicus) %>%
  write_csv('data/iSTOP/Rnorvegicus-rn6.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Rnorvegicus-rn6.csv')

# Plant ~ 400 MB
read_csv('data/CDS/Athaliana-plantsmart28.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Athaliana.TAIR.TAIR9::Athaliana, cores = cores) %>%
  locate_PAM(BSgenome.Athaliana.TAIR.TAIR9::Athaliana) %>%
  write_csv('data/iSTOP/Athaliana-plantsmart28.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Athaliana-plantsmart28.csv')

# Yeast
read_csv('data/CDS/Scerevisiae-sacCer3.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae, cores = cores) %>%
  locate_PAM(BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae) %>%
  write_csv('data/iSTOP/Scerevisiae-sacCer3.csv') %>%
  compact_iSTOP() %>%
  # Loss of cut site columns named RFLP_C_<width>
  # add_RFLP(width = 150, cores = cores) %>%
  # add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 50, cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Scerevisiae-sacCer3.csv')
