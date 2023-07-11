### Bed format SV database files ###
TRS SV database mapping by breakpoint distance (1000bp) / nonTRS SV database mapping by overlap ratio (90%) 
###########################################################################################################

1000 genome database SVs  <br />
├── 1000_gall.bp.trs                              ------1000 genomes database all trs SVs with CI info  <br />
├── 1000_gall.nontrs.gz                           ------1000 genomes database all nontrs SVs  <br /> 
├── 1000_gall.nontrs.gz.tbi    <br />
├── 1000_gall.trs                                 ------1000 genomes database all trs SVs                    
├── 1000_g.bp.trs                                 ------1000 genomes database trs SVs >=2 hits with CI info  <br />
├── 1000_g.nontrs.gz                              ------1000 genomes database nontrs SVs >=2 hits  <br />
├── 1000_g.nontrs.gz.tbi  <br />
├── 1000_g.trs                                    ------1000 genomes database trs SVs >=2 hits  <br />
├──    <br />
centromere genomic range SVs  <br />
├── centromere_qc.nontrs.gz                       ------centromere genomic range for nontrs SVs  <br />
├── centromere_qc.nontrs.gz.tbi  <br />
├── centromere_qc.trs                             ------centromere genomic range for trs SVs  <br />
├──    <br />
├── ~~(DB files normal donor WGS are considered as sequencing technical controls and not provided here)~~  <br />
├── ~~control_gall.nontrs.gz~~                        ------~~all nontrs SVs >=10 hits in normal donor WGS data~~  <br />
├── ~~control_g.nontrs.gz~~                           ------~~trs SVs >=100 hits in normal donor WGS data~~  <br />
├── ~~control_gall.trs~~                              ------~~all trs SVs >=2 hits in normal donor WGS data~~  <br />
├── ~~control_g.trs~~                                 ------~~trs SVs >=20 hits in normal donor WGS data~~  <br />
├──     <br />
COSMIC database SVs  <br />
├── cosmic_s.bp.trs                               ------COSMIC database trs SVs with CI info  <br />
├── cosmic_s.nontrs.gz                            ------COSMIC database nontrs SVs with CI info  <br />
├── cosmic_s.nontrs.gz.tbi  <br />
├── cosmic_s.trs                                   ------COSMIC database trs SVs  <br />
├──    <br />
CytoAtlas database SVs  <br />
├── cytoatlas_s.bp.trs                             ------CytoAtlas database trs SVs with CI info  <br />
├── cytoatlas_s.nontrs.gz                          ------CytoAtlas database nontrs SVs  <br />
├── cytoatlas_s.nontrs.gz.tbi  <br />
├── cytoatlas_s.trs                                ------CytoAtlas database trs SVs  <br />
├── CYTO-SV-Auto-ML.py  <br />
├── CYTO-SV-Auto-ML_tuning.py  <br />
├──    <br />
DGV database SVs  <br />
├── dgv_g.nontrs.gz                                ------DGV database nontrs SVs  <br /> 
├── dgv_g.nontrs.gz.tbi  <br />
├──    <br />
gnomAD database SVs  <br />
├── gnomad_g2ab.nontrs.gz                          ------gnomAD database nontrs SVs >=2 hits  <br />
├── gnomad_g2ab.nontrs.gz.tbi  <br />
├── gnomad_gall.nontrs.gz                          ------gnomAD database all nontrs SVs   <br />
├── gnomad_gall.nontrs.gz.tbi  <br />
├── gnomad_g.nontrs.gz                             ------gnomAD database nontrs SVs AF >=0.05  <br />
├── gnomad_g.nontrs.gz.tbi  <br />
├── gnomad_qc.nontrs.gz                            ------gnomAD database nontrs SVs with QC fail labels <br />
├── gnomad_qc.nontrs.gz.tbi  <br />
├── gnomad_qc.trs                                  ------gnomAD database trs SVs with QC fail  <br />
└── gnomad.trs                                     ------gnomAD database all trs SVs  <br />
