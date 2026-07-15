



######################################################################################################################################


# Script Remove (Filter) Transcription Factors from a vector


######################################################################################################################################
## Choose the specific set of genes


## AxonContras




##Axon Contras CPE
genesACCPE <- c("2610507B11Rik", "4632404H12Rik", "4932438A13Rik", "Accsl", "Ackr3", "Acvr1", "Ada", "Adgrv1", "Ago2", "Ago3", "Ai182371", "Aim2", "Ammecr1", "Ankrd13c", "Ankrd6", "Ano6", "Ano8", "Apln", "Arhgef1", "Arid4b", "Arl4c", "Armcx4", "Ascl1", "Atp10a", "Atp11c", "Atp2c2", "Atp8b2", "Atxn2l", "Atxn7", "B630019K06Rik", "Bahcc1", "Bcl11a", "Bcl11b", "Bcor", "Bdp1", "Blcap", "Bmpr1a", "Brwd1", "Btg1", "C2Cd3", "Cacna1b", "Cacna2d1", "Cacng7", "Cap2", "Cbln1", "Ccdc39", "Ccnl1", "Ccnt2", "Ccr9", "Cd68", "Cdc14a", "Celsr3", "Cenpi", "Cep295", "Clca2", "Clcc1", "Clcn5", "Cnot3", "Cntrl", "Creb5", "Crebzf", "Cspp1", "Ctf1", "Ctns", "Dbx1", "Dgcr8", "Diaph1", "Diaph2", "Dlg4", "Dll1", "Dlx2", "Dmxl1", "Dnah6", "Dock6", "Dyrk1A", "Dzip3", "E130311K13Rik", "E230029C05Rik", "Edem3", "Ehmt1", "Eif4g3", "Ep300", "Epc1", "Ercc5", "Erlec1", "Erp29", "Evc", "Fam13c", "Fam53a", "Fbxl5", "Fbxo11", "Fcho2", "Fgfr2", "Fktn", "Fmnl1", "Fndc3c1", "Foxd1", "Fubp3", "Gad2", "Gadd45b", "Ggta1", "Gins2", "Gja6", "Glis3", "Glrb", "Gm43517", "Gm8237", "Gpcpd1", "Gse1", "Gxylt1", "H2-t10", "Hdac4", "Hipk2", "Hs3St3b1", "Iffo1", "Ifi207", "Ikbip", "Il18bp", "Il6ra", "Inhbb", "Ino80", "Ints8", "Itga4", "Jak2", "Kcmf1", "Kcnh3", "Kcnh7", "Kdm6a", "Kdm7a", "Kmt2a", "Kmt2b", "Kmt2d", "Kmt2e", "Krba1", "Krit1", "Lipg", "Lipm", "Lix1", "Lrig2", "Lrp2", "Luc7L3", "Magel2", "Map4k2", "Map4k3", "Mbd5", "Mdn1", "Me1", "Mertk", "Mfsd4b4", "Mkrn3", "Mllt3", "Mme", "Mmp16", "Mpzl1", "Myo10", "Myo9b", "Nbeal1", "Ncor1", "Neurod2", "Nfat5", "Nfxl1", "Nin", "Nkx2-2", "Nova1", "Nox4", "Nr4a2", "Nrarp", "Nrd1", "Odf2l", "Otud4", "Pak3", "Palmd", "Pard3b", "Parp9", "Paxip1", "Pcdh10", "Pcdhb13", "Pcdhga12", "Pcnx", "Phf12", "Phf21b", "Pid1", "Pik3c2b", "Pik3cb", "Plcxd2", "Pou2f1", "Pou3f3", "Pou6f1", "Ppp4r4", "Prpf40a", "Prr14l", "Prrc2b", "Ptpn5", "Ptprz1", "Qk", "R3hdm1", "Rab23", "Rab28", "Rab31", "Rai1", "Ranbp17", "Rbbp6", "Rbm33", "Rgs7bp", "Rock2", "Rslcan18", "Rttn", "Ryk", "Samd3", "Scn1a", "Sec63", "Setd5", "Sfmbt1", "Sh2b3", "Shtn1", "Simc1", "Slc16a11", "Slc44a1", "Smg1", "Smurf2", "Sncaip", "Snrnp70", "Sntg2", "Sobp", "Sox5", "Spen", "Spopl", "Srcin1", "St18", "Stox2", "Stpg2", "Taf4", "Taok1", "Tcerg1l", "Tcf7l1", "Tchh", "Tdrp", "Tfpi", "Tmeff2", "Tmem161b", "Tmem243", "Tmem255a", "Tnip2", "Tnrc18", "Tnrc6a", "Tnrc6b", "Tnrc6c", "Tox3", "Trak1", "Trim11", "Trim59", "Trp53bp1", "Trpm7", "Ttc14", "Txlnb", "Ube2d2b", "Usp1", "Usp21", "Usp34", "Veph1", "Vezt", "Vstm2a", "Vstm2b", "Wipf3", "Zbtb18", "Zbtb25", "Zc3h12c", "Zcchc7", "Zdhhc15", "Zdhhc20", "Zfhx4", "Zfp113", "Zfp236", "Zfp280c", "Zfp280d", "Zfp334", "Zfp384", "Zfp462", "Zfp512b", "Zfp523", "Zfp609", "Zfp710")

## Axon Ipsis CPE

genesAICPE <- c("2610301B20Rik", "A430033K04Rik", "A4Galt", "Abat", "Abhd4", "Ablim3", "Acvr2a", "Adam12", "Adamtsl3", "Adcy5", "Add3", "Adgra2", "Adgrg6", "Adgrl3", "Afdn", "Aff2", "Aff3", "Aff4", "Alcam", "Ankfy1", "Ankrd12", "Ankrd50", "Anks1b", "Anxa1", "Apobec3", "Araf", "Arap2", "Arhgap12", "Arhgef12", "Arhgef6", "Armcx1", "Armcx3", "Arsb", "Arxes1", "Asic4", "Ate1", "Atp13a4", "Atp1b1", "Atp2b2", "Atp6v1a", "Atp8a1", "Atrn", "Atxn1", "Avpr1a", "Aw551984", "B4galt7", "Bend6", "Brinp3", "Brms1L", "Cacna1g", "Cacnb4", "Cadm1", "Calb1", "Calcr", "Calcrl", "Cald1", "Car10", "Carmil1", "Cask", "Cav1", "Cav2", "Cd1D1", "Cd200r1", "Cdc42ep1", "Cdk5r1", "Cers6", "Cflar", "Chic1", "Chl1", "Chn2", "Chst14", "Clcn4", "Cldnd1", "Clvs2", "Cmtm3", "Cngb1", "Colec12", "Cpne2", "Crat", "Creb3l2", "Cryzl1", "Csmd3", "Csrnp3", "Cthrc1", "Ctnnb1", "Ctsb", "Cx3cr1", "Cxadr", "Cxcl12", "Cybrd1", "Daam2", "Dact1", "Dapk2", "Dcaf7", "Dclk1", "Dgke", "Dip2a", "Dip2c", "Dmxl2", "Dnah5", "Dnajc18", "Dnmt3a", "Dock10", "Dock9", "Dpp10", "Dpp4", "Dpyd", "Ednra", "Ednrb", "Enpep", "Enpp1", "Epha4", "Ephb1", "Epm2aip1", "Esam", "Ets1", "F11r", "F2Rl2", "Fabp7", "Fam110b", "Fam120c", "Fam167a", "Fam20a", "Fam20b", "Fam214a", "Fbxl17", "Fbxo25", "Fbxo7", "Fech", "Fgfr1", "Fhl2", "Flt1", "Fmo1", "Fnbp1", "Fndc3b", "Foxc1", "Frrs1l", "Fry", "Fstl1", "Fstl5", "Fyn", "Fzd1", "Fzd3", "Gabbr1", "Gabra1", "Gabra3", "Gabrg1", "Gabrg3", "Gad1", "Gas6", "Gbp4", "Gbp9", "Gda", "Gfod1", "Gfpt1", "Gja1", "Gjd2", "Gm11733", "Gm11992", "Gm17018", "Gm5796", "Gng7", "Gnl3L", "Gnrh1", "Gpc3", "Gpm6a", "Gpm6b", "Gpr101", "Gpr146", "Gpr158", "Gpr176", "Gpr21", "Gpr26", "Gpr45", "Gpt2", "Gramd3", "Gramd4", "Gria2", "Grik1", "Grik4", "Grin2b", "Grm1", "Grm4", "Grm5", "Hbb-bt", "Hecw1", "Heg1", "Hepacam", "Herc2", "Hpcal4", "Hpgds", "Htr4", "Ica1l", "Ifi202b", "Igf2", "Igfbp5", "Igsf1", "Il6st", "Impact", "Inpp4b", "Iqub", "Itgbl1", "Itih5", "Itm2b", "Kbtbd7", "Kcnc1", "Kcnh1", "Kcnt1", "Kctd12", "Kctd12b", "Kif1b", "Kif9", "Klhl23", "Klhl32", "Klhl42", "L1cam", "Ldb2", "Lgalsl", "Lhx9", "Lima1", "Lmnb2", "Lonrf2", "Lox", "Lpp", "Lrig1", "Lrrc3b", "Lrrc55", "Lrrc8c", "Lrrn3", "Lrrtm3", "Lsamp", "Lynx1", "Mal2", "Manea", "Map3k15", "Map3k2", "Map4", "Mapk4", "Mapkbp1", "Mast4", "Mpped2", "Mr1", "Myo6", "N4bp2l1", "Nav3", "Ncam2", "Ndnf", "Neurl1b", "Nfe2L1", "Nfia", "Nfib", "Nhsl2", "Nkap", "Notch2", "Npepps", "Nptn", "Nptx2", "Npy2r", "Nr2c2", "Nrcam", "Nrp2", "Ntm", "Ntrk2", "Nwd2", "Nxph1", "Nxt2", "Oat", "Olfml2a", "Opcml", "Pabpc5", "Pcdh11x", "Pcdha12", "Pcdhac1", "Pcdhac2", "Pcdhb7", "Pcdhgc3", "Pcdhgc4", "Pcsk5", "Pde4b", "Pde7b", "Pde8b", "Pdgfrb", "Pdha1", "Pdlim5", "Pdp1", "Pex26", "Pfkfb2", "Pgap1", "Phc3", "Pkdcc", "Plat", "Plcb1", "Plcz1", "Pld2", "Plk2", "Plscr2", "Plxna4", "Plxnc1", "Ppm1f", "Ppp1r12b", "Ppp1r9a", "Prkca", "Prlr", "Prrx1", "Ptdss1", "Ptgfrn", "Pxdc1", "Pxk", "Rab11fip2", "Rab3c", "Ralyl", "Ranbp3l", "Rapgef2", "Rapgef4", "Rasa3", "Rasgrp1", "Rassf2", "Rbm20", "Reep3", "Rftn1", "Rgl1", "Rgs20", "Rgs5", "Rit2", "Rnf112", "Rnf182", "Rsf1", "S100pbp", "Samd4", "Sash1", "Scg2", "Scml4", "Selenop", "Sema5a", "Senp7", "Sez6", "Sfmbt2", "Sft2d2", "Sfxn5", "Shroom2", "Slc12a5", "Slc16a12", "Slc1a2", "Slc1a3", "Slc25a12", "Slc25a34", "Slc30a4", "Slc35f1", "Slc4a10", "Slc4a4", "Slc4a8", "Slc5a3", "Slc6a17", "Slco1a1", "Slit3", "Smad5", "Smim17", "Smurf1", "Snapc1", "Snx13", "Sorcs3", "Spag6l", "Spock3", "St8Sia2", "Stard13", "Stc1", "Strip1", "Strip2", "Syne1", "Syt1", "Syt4", "Tbc1d24", "Tbc1d2b", "Tbx15", "Tbx18", "Tcaf1", "Tcaim", "Tcf7l2", "Tenm1", "Tenm2", "Tet1", "Tex2", "Tgfb2", "Tgfbr2", "Thnsl1", "Tln2", "Tmc7", "Tmcc3", "Tmco3", "Tmem100", "Tmem106b", "Tmod2", "Tmod3", "Tmx4", "Tnik", "Tnks2", "Tns3", "Top2b", "Tpcn1", "Tpd52l1", "Tppp", "Trank1", "Tril", "Trpc4", "Tshz2", "Txnip", "Ubxn2b", "Unc5d", "Usp22", "Vangl2", "Vcam1", "Vwa3a", "Wdr17", "Wdr37", "Wdr47", "Wif1", "Wwtr1", "Xpr1", "Xrn1", "Zbtb37", "Zbtb43", "Zc3h13", "Zcchc18", "Zfp174", "Zfp275", "Zfp607b", "Zfp715", "Zfp763", "Zfp799", "Zfp941", "Zfyve16", "Zic4", "Zyg11b")

lncAC <- c(
"Xist","C130071c03rik","Gm26871","Miat","B830012l14rik","Gm26953","Gm11342",
"Gm26761","Gm26981","Gm29099","1700023c21rik","Mir124a-1hg","1700041m05rik",
"9330175e14rik","Mirg","Gm27032","A430110l20rik","Gm26629","Gm10544",
"Gm14344","Mir9-3hg","Meg3","Gm4349","Gm14664","Gm26811","Gm2694","Gm26672",
"6330403l08rik","Gm28258","Au022754","2700046a07rik","C230035i16rik",
"Gm14342","Gm15834","Gm2115","Snhg18","Gm16863","1700003m07rik",
"2610037d02rik","Gm9885","Gm28652","1700110c19rik","Gm11827","Gm10603",
"Gm12688","Thoc2l","Gm26803","Gm16892","Ino80dos","E130006d01rik","Rian",
"Mir99ahg","Firre","Smim43/gm11549","4930570d08rik","Gm15999","Snhg5",
"0610009e02rik","Trerf1","Dleu2","2610307p16rik","4930513n10rik",
"9530046b11rik","Chaserr","Snhg14","4930570g19rik","2610035f20rik",
"9430041j12rik","Tbx3os2","Gm15342","4933407k13rik","4932443l11rik",
"Gm12940","4933421o10rik","9530051g07rik","E130307a14rik","Gm20407",
"9830144p21rik","Gm4117","Gm10421","2500004c02rik","1700016p03rik",
"4632404h12rik","6530409c15rik","9330151l19rik","C330011m18rik",
"E230029c05rik","D030055h07rik","Gm15775","Gm37912","4632427e13rik",
"Gm38287","Kcnq1ot1","Kdm6bos","Gm16068","Gm13563","Sox6os","4930599n23rik",
"Gm15859",NA,"Gm14540","Bc023719","Gm14097","Gm15893","Gm37753","Dlx1as",
"Dnah2os","Gm16110","Gm15629","Itpr3os","A830011k09rik","Gm16105",
"B230398e01rik","Gm16185","Sardhos","Gm16229","Sox5os1","Gm29508",
"Gm28513","Gm15624","Gm15587","Gm15122","Gm12786","Gm4793","Gm16152",
"1700027h10rik","Gm13595","Gm12426","Gm12064","Gm17168","Gm15582",
"Uckl1os","Gm16350","A830082k12rik","Gm15813","Gad1os","Gm15423",
"A230107n01rik","4930524o07rik","9430060i03rik","Gm11670","A230077h06rik",
"Gm12531","Gm15475","Zfp85os","Gm14091","Gm17224","Gm11802","Gm11399",
"4732440d04rik","Srrm4os","Gm43331","Sox2ot","Gm42133","C030018k13rik"
)

lncAI <-  c(
"Rpph1","Rprl3","Carlr","0610040b10rik","1700007j10rik","2610300m13rik",
"2610316d01rik","2810032g03rik","3930402g23rik","4930429h19rik",
"4930553p18rik","4931403e22rik","4931406g06rik","4933431e20rik",
"6430503k07rik","6530411m01rik","9130410c08rik","9430021m05rik",
"A230020j21rik","A430027h14rik","A930012l18rik","Ai504432",
"Arhgap20os","B230216n24rik","B230334c09rik","Bc037039","Bvht",
"C030014i23rik","C230037l18rik","Celrr","D030047h15rik","D730045b01rik",
"Emx2os","Foxd2os","Gm10710","Gm11615","Gm11646","Gm12371","Gm12532",
"Gm12976","Gm13619","Gm14004","Gm14164","Gm14199","Gm14204","Gm14329",
"Gm14662","Gm15441","Gm15569","Gm15631","Gm15728","Gm15751","Gm15852",
"Gm16054","Gm16308","Gm16549","Gm17180","Gm17322","Gm19590","Gm20412",
"Gm20528","Gm20559","Gm2415","Gm26691","Gm26716","Gm26813","Gm26982",
"Gm26995","Gm28289","Gm28892","Gm29587","Gm29681","Gm33051","Gm37233",
"Gm37756","Gm38165","Gm38336","Gm42488","Gm42808","Gm42836","Gm44317",
"Gm4675","Gm9831","Gm9866","H19","Halr1","Lhx1os","Lncenc1","Mexis",
"Pcsk2os2","R74862","Rapgef4os2","Zim3","4930533k18rik","Gm10308","Gm10382"
)
)

## Choose your option> 

# genes <- genesACCPE
# genes <- genesAICPE
# genes <- lncAC
genes <- lncAI

######################################################################################################################################

library(org.Mm.eg.db)
library(AnnotationDbi)
library(tidyverse)

tf_entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = keys(org.Mm.eg.db, keytype = "GOALL"),
  columns = c("SYMBOL"),
  keytype = "GOALL"
) %>%
  filter(GOALL == "GO:0003700") %>%   # transcription factor activity
  pull(SYMBOL) %>%
  unique()
  
  
genes_tf <- intersect(genes, tf_entrez)
genes_no_tf <- setdiff(genes, tf_entrez)


length(genes)        # original
length(genes_tf)     # cuántos TF quitaste
length(genes_no_tf)  # limpio

## Axon contras

# [1] 282
# [1] 35
# [1] 247

## Axon Ipsis

# [1] 433
# [1] 27
# [1] "A430033K04Rik" "Aff3"          "Creb3l2"       "Csrnp3"       
# [5] "Ets1"          "Foxc1"         "Lhx9"          "Nfia"         
# [9] "Nfib"          "Nr2c2"         "Plscr2"        "Prrx1"        
# [13] "Smad5"         "Tbx15"         "Tbx18"         "Tcf7l2"       
# [17] "Tshz2"         "Zbtb37"        "Zbtb43"        "Zfp174"       
# [21] "Zfp275"        "Zfp607b"       "Zfp715"        "Zfp763"       
# [25] "Zfp799"        "Zfp941"        "Zic4"      

# [1] 406


###################################################################################################

library(clusterProfiler)
library(ReactomePA)
library(enrichR)
# library(org.Hs.eg.db)   # cambia si es mouse
library(org.Mm.eg.db)   # cambia si es mouse
library(tidyverse)

## Convertir a Entrez (clave para clusterProfiler)

gene_df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

entrez <- gene_df$ENTREZID
 # Axon Contras 3.9% of input gene IDs are fail to map...
 # Axon Ipsis  2.77% of input gene IDs are fail to map...


gene_dfnoTF <- bitr(
  genes_no_tf,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Mm.eg.db
)

entrez2 <- gene_dfnoTF$ENTREZID
 # Axon Contras 4.45% of input gene IDs are fail to map...
 # Axon Ipsis  2.96% of input gene IDs are fail to map...

#############################################################

### CPE 

# GO Biological Process (clusterProfiler)
ego_bp <- enrichGO(
  gene          = entrez,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)

dotplot(ego_bp, showCategory = 20)
ggsave("GOBP_CPEAxons.svg")



library(enrichplot)
library(ggplot2)

p <- dotplot(ego_bp, showCategory = 20)

# p + 
 #  scale_color_gradient(
  #   low = "yellow",
   # high = "red"
  # )

p + 
  scale_color_gradientn(
    colours = c("yellow", "orange", "red")
  )


ggsave("GOBP_CPEAxons.svg", width = 8, height = 6)  

res <- as.data.frame(ego_bp)

res_top <- res %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(
    logFDR = -log10(p.adjust),
    Pathway = stringr::str_wrap(Description, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
  
  ggplot(res_top,
       aes(x = logFDR,
           y = Pathway)) +

  geom_segment(aes(x = 0,
                   xend = logFDR,
                   yend = Pathway),
               color = "grey60",
               linewidth = 0.9) +

  geom_point(aes(size = Count,
                 color = logFDR)) +

  scale_color_gradient(low = "yellow",
                       high = "red") +

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene count",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  ggsave("GOBP_CPEAxons_lolliplot.svg")

## CPE noTF
ego_bp2 <- enrichGO(
  gene          = entrez2,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  readable      = TRUE
)


dotplot(ego_bp2, showCategory = 20)
ggsave("GOBP_CPEAxons_noTF.svg")


res2 <- as.data.frame(ego_bp2)

res_top2 <- res2 %>%
  arrange(p.adjust) %>%
  slice_head(n = 15) %>%
  mutate(
    logFDR = -log10(p.adjust),
    Pathway = stringr::str_wrap(Description, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
  
  ggplot(res_top2,
       aes(x = logFDR,
           y = Pathway)) +

  geom_segment(aes(x = 0,
                   xend = logFDR,
                   yend = Pathway),
               color = "grey60",
               linewidth = 0.9) +

  geom_point(aes(size = Count,
                 color = logFDR)) +

  scale_color_gradient(low = "yellow",
                       high = "red") +

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene count",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  ggsave("GOBP_CPEAxonsnoTF_lolliplot.svg")
  
###############################################################

## Reactome 

library(ReactomePA)
react <- enrichPathway(
  gene          = entrez,
  organism      = "mouse",
  pAdjustMethod = "BH",
  readable      = TRUE
)

dotplot(react, showCategory = 20)
ggsave("reactome_CPEAxons.svg")


react
#
# over-representation test
#
#...@organism 	 mouse 
#...@ontology 	 Reactome 
#...@keytype 	 ENTREZID 
#...@gene 	 chr [1:271] "74034" "381411" "12778" "11477" "11486" "110789" "239528" ...
#...pvalues adjusted by 'BH' with cutoff <0.05 
#...0 enriched terms found
#...Citation
  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
  reactome pathway analysis and visualization. Molecular BioSystems
  2016, 12(2):477-479 


 

react2 <- enrichPathway(
  gene          = entrez2,
  organism      = "mouse",
  pAdjustMethod = "BH",
  readable      = TRUE
)
 
dotplot(react2, showCategory = 20)
ggsave("reactome_CPEAxons_noTF.svg")

over-representation test
#
#...@organism 	 mouse 
#...@ontology 	 Reactome 
#...@keytype 	 ENTREZID 
#...@gene 	 chr [1:236] "74034" "381411" "12778" "11477" "11486" "110789" "239528" ...
#...pvalues adjusted by 'BH' with cutoff <0.05 
#...0 enriched terms found
#...Citation
  Guangchuang Yu, Qing-Yu He. ReactomePA: an R/Bioconductor package for
  reactome pathway analysis and visualization. Molecular BioSystems
  2016, 12(2):477-479 

######################################################################################################################################

library(enrichR)
Welcome to enrichR
Checking connections ... 
Enrichr ... Connection is Live!
FlyEnrichr ... Connection is Live!
WormEnrichr ... Connection is Live!
YeastEnrichr ... Connection is Live!
FishEnrichr ... Connection is Live!
OxEnrichr ... Connection is Live!
> dbs <- enrichR::listEnrichrDbs()
> dbs
 
                                           libraryName
1                                  Genome_Browser_PWMs
2                             TRANSFAC_and_JASPAR_PWMs
3                            Transcription_Factor_PPIs
4                                            ChEA_2013
5                     Drug_Perturbations_from_GEO_2014
6                              ENCODE_TF_ChIP-seq_2014
7                                        BioCarta_2013
8                                    WikiPathways_2013
9                  Disease_Signatures_from_GEO_up_2014
10                                           KEGG_2013
11                          TF-LOF_Expression_from_GEO
12                                 TargetScan_microRNA
13                                    PPI_Hub_Proteins
14                                           GeneSigDB
15                                 Chromosome_Location
16                                    Human_Gene_Atlas
17                                    Mouse_Gene_Atlas
18                            Human_Phenotype_Ontology
19                     Epigenomics_Roadmap_HM_ChIP-seq
20                                            KEA_2013
21                   NURSA_Human_Endogenous_Complexome
22                                               CORUM
23                             SILAC_Phosphoproteomics
24                                         Old_CMAP_up
25                                       Old_CMAP_down
26                                        OMIM_Disease
27                                       OMIM_Expanded
28                                           VirusMINT
29                                MSigDB_Computational
30                         MSigDB_Oncogenic_Signatures
31               Disease_Signatures_from_GEO_down_2014
32                     Virus_Perturbations_from_GEO_up
33                   Virus_Perturbations_from_GEO_down
34                       Cancer_Cell_Line_Encyclopedia
35                            NCI-60_Cancer_Cell_Lines
36         Tissue_Protein_Expression_from_ProteomicsDB
37   Tissue_Protein_Expression_from_Human_Proteome_Map
38                                    HMDB_Metabolites
39                               Pfam_InterPro_Domains
40                                Allen_Brain_Atlas_up
41                             ENCODE_TF_ChIP-seq_2015
42                   ENCODE_Histone_Modifications_2015
43                   Phosphatase_Substrates_from_DEPOD
44                              Allen_Brain_Atlas_down
45                   ENCODE_Histone_Modifications_2013
46                           Achilles_fitness_increase
47                           Achilles_fitness_decrease
48                                       BioCarta_2015
49                                       HumanCyc_2015
50                                           KEGG_2015
51                                        Panther_2015
52                                   WikiPathways_2015
53                                              ESCAPE
54                                          HomoloGene
55                 Disease_Perturbations_from_GEO_down
56                   Disease_Perturbations_from_GEO_up
57                    Drug_Perturbations_from_GEO_down
58                    Genes_Associated_with_NIH_Grants
59                      Drug_Perturbations_from_GEO_up
60                                            KEA_2015
61                      Gene_Perturbations_from_GEO_up
62                    Gene_Perturbations_from_GEO_down
63                                           ChEA_2015
64                                               dbGaP
65                            LINCS_L1000_Chem_Pert_up
66                          LINCS_L1000_Chem_Pert_down
67                         GTEx_Tissue_Expression_Down
68                           GTEx_Tissue_Expression_Up
69                  Ligand_Perturbations_from_GEO_down
70                   Aging_Perturbations_from_GEO_down
71                     Aging_Perturbations_from_GEO_up
72                    Ligand_Perturbations_from_GEO_up
73                    MCF7_Perturbations_from_GEO_down
74                      MCF7_Perturbations_from_GEO_up
75                 Microbe_Perturbations_from_GEO_down
76                   Microbe_Perturbations_from_GEO_up
77               LINCS_L1000_Ligand_Perturbations_down
78                 LINCS_L1000_Ligand_Perturbations_up
79            L1000_Kinase_and_GPCR_Perturbations_down
80              L1000_Kinase_and_GPCR_Perturbations_up
81                                           KEGG_2016
82                                   WikiPathways_2016
83           ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X
84                  Kinase_Perturbations_from_GEO_down
85                    Kinase_Perturbations_from_GEO_up
86                                       BioCarta_2016
87                                       HumanCyc_2016
88                                     NCI-Nature_2016
89                                        Panther_2016
90                                          DrugMatrix
91                                           ChEA_2016
92                                               huMAP
93                                      Jensen_TISSUES
94   RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO
95                                 Jensen_COMPARTMENTS
96                                     Jensen_DISEASES
97                                        BioPlex_2017
98                                      ARCHS4_Tissues
99                                   ARCHS4_Cell-lines
100                                   ARCHS4_IDG_Coexp
101                               ARCHS4_Kinases_Coexp
102                                   ARCHS4_TFs_Coexp
103                            SysMyo_Muscle_Gene_Sets
104                                    miRTarBase_2017
105                           TargetScan_microRNA_2017
106               Enrichr_Libraries_Most_Popular_Genes
107            Enrichr_Submissions_TF-Gene_Coocurrence
108         Data_Acquisition_Method_Most_Popular_Genes
109                                             DSigDB
110            TF_Perturbations_Followed_by_Expression
111                           Chromosome_Location_hg19
112           Rare_Diseases_AutoRIF_ARCHS4_Predictions
113           Rare_Diseases_GeneRIF_ARCHS4_Predictions
114                   Rare_Diseases_GeneRIF_Gene_Lists
115                   Rare_Diseases_AutoRIF_Gene_Lists
116                                    SubCell_BarCode
117                                  GWAS_Catalog_2019
118                            WikiPathways_2019_Human
119                            WikiPathways_2019_Mouse
120                  TRRUST_Transcription_Factors_2019
121                                    KEGG_2019_Human
122                                    KEGG_2019_Mouse
123                              InterPro_Domains_2019
124                                  Pfam_Domains_2019
125      DepMap_WG_CRISPR_Screens_Broad_CellLines_2019
126     DepMap_WG_CRISPR_Screens_Sanger_CellLines_2019
127                                 UK_Biobank_GWAS_v1
128                                     BioPlanet_2019
129                                       ClinVar_2019
130                                        PheWeb_2019
131                                           DisGeNET
132                               HMS_LINCS_KinomeScan
133                               CCLE_Proteomics_2020
134                                  ProteomicsDB_2020
135                        lncHUB_lncRNA_Co-Expression
136                      Virus-Host_PPI_P-HIPSTer_2020
137                        Elsevier_Pathway_Collection
138                     Table_Mining_of_CRISPR_Studies
139                         COVID-19_Related_Gene_Sets
140                               MSigDB_Hallmark_2020
141               Enrichr_Users_Contributed_Lists_2020
142                                      TG_GATES_2020
143                   Allen_Brain_Atlas_10x_scRNA_2021
144               Descartes_Cell_Types_and_Tissue_2021
145                                    KEGG_2021_Human
146                             WikiPathway_2021_Human
147 HuBMAP_ASCT_plus_B_augmented_w_RNAseq_Coexpression
148                         GO_Biological_Process_2021
149                         GO_Cellular_Component_2021
150                         GO_Molecular_Function_2021
151               MGI_Mammalian_Phenotype_Level_4_2021
152                          CellMarker_Augmented_2021
153                            Orphanet_Augmented_2021
154                    COVID-19_Related_Gene_Sets_2021
155                           PanglaoDB_Augmented_2021
156                            Azimuth_Cell_Types_2021
157                          PhenGenI_Association_2021
158                         GTEx_Aging_Signatures_2021
159                                 HDSigDB_Human_2021
160                                 HDSigDB_Mouse_2021
161                    HuBMAP_ASCTplusB_augmented_2022
162                             FANTOM6_lncRNA_KD_DEGs
163                           MAGMA_Drugs_and_Diseases
164                                     Tabula_Sapiens
165                                          ChEA_2022
166                    Diabetes_Perturbations_GEO_2022
167               LINCS_L1000_Chem_Pert_Consensus_Sigs
168               LINCS_L1000_CRISPR_KO_Consensus_Sigs
169                                       Tabula_Muris
170                                      Reactome_2022
171                                         SynGO_2022
172                  GlyGen_Glycosylated_Proteins_2022
173                              IDG_Drug_Targets_2022
174                        KOMP2_Mouse_Phenotypes_2022
175            Metabolomics_Workbench_Metabolites_2022
176                         Proteomics_Drug_Atlas_2023
177                            The_Kinase_Library_2023
178                               GTEx_Tissues_V8_2023
179                         GO_Biological_Process_2023
180                         GO_Cellular_Component_2023
181                         GO_Molecular_Function_2023
182                                PFOCR_Pathways_2023
183                                  GWAS_Catalog_2023
184                                      GeDiPNet_2023
185                                        MAGNET_2023
186                                       Azimuth_2023
187                                  Rummagene_kinases
188                               Rummagene_signatures
189                    Rummagene_transcription_factors
190                                       MoTrPAC_2023
191                             WikiPathway_2023_Human
192        DepMap_CRISPR_GeneDependency_CellLines_2023
193                                         SynGO_2024
194                                    CellMarker_2024
195                            DGIdb_Drug_Targets_2024
196               MGI_Mammalian_Phenotype_Level_4_2024
197                            WikiPathways_2024_Mouse
198                            WikiPathways_2024_Human
199                            The_Kinase_Library_2024
200                                       PerturbAtlas
201                             Reactome_Pathways_2024
202                         GO_Biological_Process_2025
203                         GO_Cellular_Component_2025
204                         GO_Molecular_Function_2025
205                               NIBR_DRUGseq_2025_up
206                             NIBR_DRUGseq_2025_down
207                       Jensen_DISEASES_Curated_2025
208                  Jensen_DISEASES_Experimental_2025
209                          COMPARTMENTS_Curated_2025
210             PerturbAtlas_MouseGenePerturbationSigs
211                            PerturbSeq_ReplogleK562
212                          TISSUES_Experimental_2025
213                     COMPARTMENTS_Experimental_2025
214                            PerturbSeq_ReplogleRPE1
215                               TISSUES_Curated_2025
216                    RummaGEO_DrugPerturbations_2025
217                    RummaGEO_GenePerturbations_2025
218                                  GWAS_Catalog_2025
219                                       ClinVar_2025
220         CM4AI_U2OS_Protein_Localization_Assemblies
221                              JASPAR_PWM_Human_2025
222                              JASPAR_PWM_Mouse_2025
223          Sciplex_Drug_Perturbation_Signatures_2025
224                                      Carcinogenome
225                                          KEGG_2026


dbs_use <- c(
  "WikiPathways_2019_Mouse",
"WikiPathways_2024_Mouse",
"KEGG_2019_Mouse",
"JASPAR_PWM_Mouse_2025",
"HDSigDB_Mouse_2021",
"TargetScan_microRNA",
"Reactome_Pathways_2024",
"GO_Biological_Process_2023",
"lncHUB_lncRNA_Co-Expression",
"Virus-Host_PPI_P-HIPSTer_2020",
"Elsevier_Pathway_Collection",
  "NCI-Nature_2016",
  "Panther_2016"
)

enrichr_res <- enrichr(genes, dbs_use)
enrichr_resnoTF <- enrichr(genes_no_tf, dbs_use)



#######################################################################################

## Panther_2016
       

Panther_2016 <- enrichr_res[["Panther_2016"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- Panther_2016 %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("Panther_2016.svg")
       write.table(Panther_2016, "Panther_2016.txt" )
  
  
 Panther_2016noTF <- enrichr_resnoTF[["Panther_2016"]] %>%
  arrange(Adjusted.P.value)



wp_top <- Panther_2016noTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("Panther_2016noTF.svg")
         write.table(Panther_2016noTF, "Panther_2016noTF.txt" )
         
         
#######################################################################################

## Elsevier_Pathway_Collection
Elsevier_Pathway_Collection <- enrichr_res[["Elsevier_Pathway_Collection"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- Elsevier_Pathway_Collection %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("Elsevier_Pathway_Collection.svg")
       write.table(Elsevier_Pathway_Collection, "Elsevier_Pathway_Collection.txt" )
  
  
 Elsevier_Pathway_CollectionnoTF <- enrichr_resnoTF[["Elsevier_Pathway_Collection"]] %>%
  arrange(Adjusted.P.value)



wp_top <- Elsevier_Pathway_CollectionnoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("Elsevier_Pathway_CollectionnoTF.svg")
         write.table(Elsevier_Pathway_CollectionnoTF, "Elsevier_Pathway_CollectionnoTF.txt" )
         
         
#######################################################################################

## NCINature_2016
NCINature_2016 <- enrichr_res[["NCI-Nature_2016"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- NCINature_2016 %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("NCINature_2016.svg")
       write.table(NCINature_2016, "NCINature_2016.txt" )
  
  
 NCINature_2016noTF <- enrichr_resnoTF[["NCI-Nature_2016"]] %>%
  arrange(Adjusted.P.value)



wp_top <- NCINature_2016noTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("NCINature_2016noTF.svg")
         write.table(NCINature_2016noTF, "NCINature_2016noTF.txt" )
         
         
#######################################################################################

## lncHUB_lncRNA_Co-Expression 
lncHUB_lncRNA_CoExpression <- enrichr_res[["lncHUB_lncRNA_Co-Expression"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- lncHUB_lncRNA_CoExpression %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("lncHUB_lncRNA_CoExpression.svg")
       write.table(lncHUB_lncRNA_CoExpression, "lncHUB_lncRNA_CoExpression.txt" )
  
  
 lncHUB_lncRNA_CoExpressionnoTF <- enrichr_resnoTF[["lncHUB_lncRNA_Co-Expression"]] %>%
  arrange(Adjusted.P.value)



wp_top <- lncHUB_lncRNA_CoExpressionnoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("lncHUB_lncRNA_CoExpressionnoTF.svg")
         write.table(lncHUB_lncRNA_CoExpressionnoTF, "lncHUB_lncRNA_CoExpressionnoTF.txt" )
         
         
#######################################################################################

## JASPAR_PWM_Mouse_2025 
JASPAR_PWM_Mouse_2025 <- enrichr_res[["JASPAR_PWM_Mouse_2025"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- JASPAR_PWM_Mouse_2025 %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("JASPAR_PWM_Mouse_2025.svg")
       write.table(JASPAR_PWM_Mouse_2025, "JASPAR_PWM_Mouse_2025.txt" )
  
  
 JASPAR_PWM_Mouse_2025noTF <- enrichr_resnoTF[["JASPAR_PWM_Mouse_2025"]] %>%
  arrange(Adjusted.P.value)



wp_top <- JASPAR_PWM_Mouse_2025noTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("JASPAR_PWM_Mouse_2025noTF.svg")
         write.table(JASPAR_PWM_Mouse_2025noTF, "JASPAR_PWM_Mouse_2025noTF.txt" )
         
#######################################################################################

## TargetScan_microRNA 
TargetScan_microRNA <- enrichr_res[["TargetScan_microRNA"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- TargetScan_microRNA %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("TargetScan_microRNA.svg")
       write.table(TargetScan_microRNA, "TargetScan_microRNA.txt" )
  
  
 TargetScan_microRNAnoTF <- enrichr_resnoTF[["TargetScan_microRNA"]] %>%
  arrange(Adjusted.P.value)



wp_top <- TargetScan_microRNAnoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("TargetScan_microRNAnoTF.svg")
         write.table(TargetScan_microRNAnoTF, "TargetScan_microRNAnoTF.txt" )
         
#######################################################################################


## GOBIOLOGICAL PROCESS
gobp <- enrichr_res[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- gobp %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("gobp_ACCPE.svg")
    
 write.table(gobp, "GOBP_ACCPEtxt" )
  
  library(tidyverse)
gobpnoTF <- enrichr_resnoTF[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value)

wp_top <- gobpnoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("gobpnoTF_ACCPE_noTF.svg")
  
 write.table(gobpnoTF, "GOBP_ACCPE_NOtf.txt" )




#######################################################################################



## KEGG_2019_Mouse
reacto <- enrichr_res[["Reactome_Pathways_2024"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- reacto %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("reactome_ACCPE.svg")
   write.table(reacto, "reactome_ACCPE.txt" )
  
  library(tidyverse)
reactonoTF <- enrichr_resnoTF[["Reactome_Pathways_2024"]] %>%
  arrange(Adjusted.P.value)

wp_top <- reactonoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("reactomenoTF_ACCPE_noTF.svg")
     write.table(reactonoTF, "reactome_ACCPE_noTF.txt" )




#######################################################################################

## KEGG_2019_Mouse
kegg <- enrichr_res[["KEGG_2019_Mouse"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- kegg %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("kegg_ACCPE.svg")
     write.table(kegg, "kegg_ACCPE.txt" )
  
  
  library(tidyverse)
 keggnoTF <- enrichr_resnoTF[["KEGG_2019_Mouse"]] %>%
  arrange(Adjusted.P.value)

wp_top <- keggnoTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("keggnoTF_ACCPE_noTF.svg")
    write.table(keggnoTF, "keggnoTF_ACCPE_noTF.txt" )

#######################################################################################

## WikiPathways 
wp2024 <- enrichr_res[["WikiPathways_2024_Mouse"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- wp2024 %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("WP2024_ACCPE.svg")
       write.table(wp2024, "WP2024_ACCPE.txt" )
  
  
 wp2024noTF <- enrichr_resnoTF[["WikiPathways_2024_Mouse"]] %>%
  arrange(Adjusted.P.value)

   
  library(tidyverse)

wp2024noTF <- enrichr_resnoTF[["WikiPathways_2024_Mouse"]] %>%
  arrange(Adjusted.P.value)

wp_top <- wp2024noTF %>%
  filter(!is.na(Adjusted.P.value)) %>%
  separate(Overlap,
           into = c("Hits", "SetSize"),
           sep = "/",
           convert = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )

ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +
  geom_point() +
  scale_color_gradient(low = "yellow", high = "red") +  # skyblue
  theme_classic(base_size = 14)

  ggsave("WP2024_ACCPE_noTF.svg")
         write.table(wp2024noTF, "WP2024_ACCPE_noTF.txt" )
         
######################################################################################################################################
## 

featureCounts -T 32 -a gencode.vM38.annotation.gtf -o counts.txt *.bam
featureCounts -T 32 -s 2 -t exon -g gene_name -a gencode.vM38.long_noncoding_RNAs.gtf -o counts_lncRNA.txt *.bam


library(dplyr)

counts <- read.delim("mRNATPMcounts.txt", header = TRUE, stringsAsFactors = FALSE)
mRNA_TPM <- counts
# Eliminar genes con expresión cero en todas las muestras
mRNA_TPM <- mRNA_TPM[rowSums(mRNA_TPM[,-1]) > 0, ]

# Poner GENEID como rownames y quedarte solo con la matriz de TPM
rownames(mRNA_TPM) <- mRNA_TPM$GENEID
mRNA_TPM <- mRNA_TPM[,-1]

# Log-transform (opcional, útil para correlaciones)
logTPM_mRNA <- log2(mRNA_TPM + 1)

lncRNA_TPM <- read.delim("lncRNAcountsTPM.txt", header = TRUE, stringsAsFactors = FALSE)

TPM_matrix <- lncRNA_TPM[, grep("_sort.bam", colnames(lncRNA_TPM))]
rownames(TPM_matrix) <- lncRNA_TPM$GENEID

logTPM <- log2(TPM_matrix + 1)  # Normalization 


cor_mat <- cor(t(logTPM), t(logTPM_mRNA), method = "spearman")
targets <- names(which(cor_mat["Xist", ] > 0.4))

# asegurarse de que las muestras estén en el mismo orden
all(colnames(logTPM) == colnames(logTPM_mRNA))  # debería dar TRUE

# Si no da TRUE, reordena:
# logTPM_mRNA <- logTPM_mRNA[, colnames(logTPM)]
# matriz de correlaciones: filas = lncRNA, columnas = mRNA
cor_mat <- cor(t(logTPM), t(logTPM_mRNA), method = "spearman")


# Extraer pares de interés
# Por ejemplo, todos los pares con correlación mayor a 0.4 o menor a -0.4:

library(tidyr)

# convertir a tabla larga
cor_long <- as.data.frame(as.table(cor_mat))
colnames(cor_long) <- c("lncRNA", "mRNA", "correlation")

# filtrar solo correlaciones fuertes
cor_pairs <- cor_long %>% filter(abs(correlation) > 0.4)

head(cor_pairs)
         lncRNA    mRNA correlation
1       Gm36635 Gm38148  -0.4803845
2       Gm16675 Gm38148  -0.4803845
3       Gm35835 Gm38148   0.4803845
4 E130307A14Rik Gm38148   0.4803845
5       Atcayos Gm38148   0.4803845
6       Gm10754 Gm38148   0.4803845

 write.csv(cor_pairs, "lncRNA_mRNA_correlations.csv", row.names = FALSE)


library(dplyr)

# cuántos mRNAs correlaciona cada lncRNA
lncRNA_summary <- cor_pairs %>%
  group_by(lncRNA) %>%
  summarise(n_targets = n(), 
            mean_cor = mean(correlation))

head(lncRNA_summary)
head(lncRNA_summary)
# A tibble: 6 × 3
  lncRNA  n_targets mean_cor
  <fct>       <int>    <dbl>
1 Gm60430      1173   -0.308
2 Gm38947      1154    0.265
3 Gm49291      1155   -0.291
4 Ttc24        1130   -0.299
5 Gm72405       784    0.122
6 Gm57403       702    0.105


mRNA_summary <- cor_pairs %>%
  group_by(mRNA) %>%
  summarise(n_lncRNA = n(), mean_cor = mean(correlation))

head(mRNA_summary)
# Sirve para ver qué genes codificantes están más ligados a lncRNAs, ideal si luego quieres hacer GO enrichment.
# A tibble: 6 × 3
  mRNA               n_lncRNA mean_cor
  <fct>                 <int>    <dbl>
1 Gm38148                2187    0.152
2 Gm38385                1912    0.278
3 Gm70753                1620    0.335
4 ENSMUSG00000131850     3263    0.196
5 Gm37144                3561    0.273
6 A930006A01Rik          2406    0.237

library(pheatmap)

# Por memoria, quizá filtrar solo top lncRNAs y genes más conectados
top_lnc <- lncRNA_summary %>% top_n(20, n_targets) %>% pull(lncRNA)
top_mRNA <- mRNA_summary %>% top_n(30, n_lncRNA) %>% pull(mRNA)

heatmap_matrix <- cor_mat[top_lnc, top_mRNA]
pheatmap(heatmap_matrix)

#############################################################################################

# Analizar diferencias entre E13.5 y E15.5

# Si quieres comparar retina vs CH en distintos tiempos, podrías crear un vector de condiciones:

samples <- colnames(logTPM)  # o logTPM_mRNA
condition <- ifelse(grepl("E13_5_CH", samples), "CH_E13",
             ifelse(grepl("E15_5_CH", samples), "CH_E15",
             ifelse(grepl("E13_5_retina", samples), "Retina_E13",
             "Retina_E15")))
             
library(dplyr)

CH13 = contralateral Axons
CH15 = Ipsilateral aXons 
# Ejemplo para lncRNAs
avg_expr <- data.frame(lncRNA = rownames(logTPM))

avg_expr$CH_E13 <- rowMeans(logTPM[, condition=="CH_E13"])
avg_expr$CH_E15 <- rowMeans(logTPM[, condition=="CH_E15"])

# Fold change CH E15 vs CH E13
avg_expr$FC_E15_E13 <- avg_expr$CH_E15 - avg_expr$CH_E13  # log2(TPM+1), así que resta = log2FC

E13_specific <- avg_expr$lncRNA[avg_expr$FC_E15_E13 < -1]  # más de 2x en E13
E15_specific <- avg_expr$lncRNA[avg_expr$FC_E15_E13 > 1]   # más de 2x en E15

# Crear factor y asignar colores
lncRNA_sidecolor <- factor(lncRNA_sidecolor, 
                           levels = c("Contralateral_Axons", "Ipsilateral_Axons", "Inespecific"))

ann_colors <- list(Specificity = c(
  Contralateral_Axons = "red",
  Ipsilateral_Axons = "skyblue",
  gray = "gray"
))

pheatmap(heatmap_matrix, 
         annotation_row = data.frame(Specificity = lncRNA_sidecolor),
         annotation_colors = ann_colors)
         
         ################################################
         
         
# filtrar correlaciones fuertes
cor_pairs_filtered <- cor_pairs %>% filter(abs(correlation) > 0.4)

# ejemplo para un lncRNA específico, o todos juntos
targets <- unique(cor_pairs_filtered$mRNA)


library(clusterProfiler)
library(org.Mm.eg.db)  # ratón

ego <- enrichGO(
  gene          = targets,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",    # si tus mRNA están en SYMBOL
  ont           = "BP",        # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# ver resultados
head(ego)
head(ego)
[1] ID          Description GeneRatio   BgRatio     pvalue      p.adjust   
[7] qvalue      geneID      Count      
<0 rows> (or 0-length row.names)


# Barplot de top 20 procesos
barplot(ego, showCategory = 20)

# O dotplot
dotplot(ego, showCategory = 20)

         library(biomaRt)

# Conectar a Ensembl Mouse
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# vector de tus genes RIK/Gm
genes_ensembl <- targets  

# mapear a SYMBOL
gene_map <- getBM(
  filters = "external_gene_name",    # si tus nombres son RIK/Gm, este filtro funciona
  attributes = c("external_gene_name", "ensembl_gene_id"),
  values = genes_ensembl,
  mart = mart
)

# ahora usa gene_map$external_gene_name en enrichGO
targets_symbol <- gene_map$external_gene_name

ego <- enrichGO(
  gene = targets_symbol,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

head(ego)
head(ego)
[1] ID          Description GeneRatio   BgRatio     pvalue      p.adjust   
[7] qvalue      geneID      Count      
<0 rows> (or 0-length row.names)

length(targets_symbol)
length(unique(targets_symbol))

[1] 3638
[1] 3618

# Obtener Ensembl IDs
gene_map <- getBM(
  filters = "external_gene_name",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = targets,
  mart = mart
)

ensembl_ids <- gene_map$ensembl_gene_id

# enrichGO usando Ensembl IDs
ego <- enrichGO(
  gene = ensembl_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

head(ego)

#############################################################################


near to lncRNA 

library(rtracklayer)
gtf <- import("gencode.vM38.annotation.gtf")

# filtrar genes
genes <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes)

genes <- gtf[gtf$type == "gene"]
genes_df <- as.data.frame(genes)

# Protein-coding y lncRNAs
protein_coding <- genes_df %>% filter(gene_type == "protein_coding")
lncRNA_all <- genes_df %>% filter(gene_type %in% c("lncRNA", "3prime_overlapping_ncRNA", "antisense", "lincRNA"))

lncAC_df <- lncRNA_all %>% filter(gene_name %in% lncAC) # Filtrar solo tus lncRNAs de interés
lncAI_df <- lncRNA_all %>% filter(gene_name %in% lncAI) # Filtrar solo tus lncRNAs de interés
gr_genes <- GRanges(seqnames = protein_coding$seqnames,
                    ranges = IRanges(start = protein_coding$start,
                                     end = protein_coding$end),
                    gene = protein_coding$gene_name)

gr_lnc <- GRanges(seqnames = lncAC_df$seqnames,
                  ranges = IRanges(start = lncAC_df$start,
                                   end = lncAC_df$end),
                  gene = lncAC_df$gene_name)       
                  
gr_lncAI <- GRanges(seqnames = lncAI_df$seqnames,
                  ranges = IRanges(start = lncAI_df$start,
                                   end = lncAI_df$end),
                  gene = lncAI_df$gene_name)  

closest_idx <- nearest(gr_lnc, gr_genes)
lncAC_df$closest_gene <- gr_genes$gene[closest_idx] 

closest_idxAI <- nearest(gr_lncAI, gr_genes)
lncAI_df$closest_gene <- gr_genes$gene[closest_idxAI] 



library(dplyr)

lncAC_df %>% dplyr::select(gene_name, closest_gene) %>% head(20) 
lncAI_df %>% dplyr::select(gene_name, closest_gene) %>% head(20) 
closest_genes <- unique(lncAC_df$closest_gene)
closest_genes <- closest_genes[!is.na(closest_genes)]  # eliminar posibles NA

closest_genesAI <- unique(lncAI_df$closest_gene)
closest_genes <- closest_genesAI[!is.na(closest_genesAI)]  # eliminar posibles NA

library(clusterProfiler)
library(org.Mm.eg.db)

ego <- enrichGO(
  gene          = closest_genes,   # tus genes codificantes
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",        # si tus genes son SYMBOLs
  ont           = "BP",            # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

# ver resultados
head(ego)    
write.table("GOBP_neargenelncRNA.txt")
# Barplot de los top 20 procesos
barplot(ego, showCategory = 20)

# Dotplot (más compacto)
dotplot(ego, showCategory = 20)
                 
                 
enrichr_res <- enrichr(closest_genes, dbs_use)

## WikiPathways 
WP_neargenelncRNA <- enrichr_res[["WikiPathways_2024_Mouse"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- WP_neargenelncRNA %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("WP_neargenelncRNA.svg")
       write.table(WP_neargenelncRNA, "WP_neargenelncRNA.txt" )
  
  
 ## Reactome_Pathways_2024 
Reactome_neargenelncRNA <- enrichr_res[["Reactome_Pathways_2024"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- Reactome_neargenelncRNA %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("Reactome_neargenelncRNA.svg")
       write.table(Reactome_neargenelncRNA, "Reactome_neargenelncRNA.txt" )
  
  
  ## GO_Biological_Process_2023 
GOBP_neargenelncRNA <- enrichr_res[["GO_Biological_Process_2023"]] %>%
  arrange(Adjusted.P.value)

library(tidyverse)

wp_top <- GOBP_neargenelncRNA %>%
  filter(!is.na(Adjusted.P.value)) %>%
  slice_head(n = 20) %>%
  mutate(
    logFDR = -log10(Adjusted.P.value),
    Pathway = stringr::str_wrap(Term, width = 40),
    Pathway = forcats::fct_reorder(Pathway, logFDR)
  )
 
  
  wp_top <- wp_top %>%
  separate(Overlap, into = c("Hits", "SetSize"), sep = "/", convert = TRUE)
  
  ggplot(wp_top,
       aes(x = logFDR,
           y = Pathway,
           size = Hits,
           color = logFDR)) +

  geom_point() +

  scale_color_gradient(low = "yellow",
                       high = "red") +     # skyblue

  labs(
    x = "-log10(adj p-value)",
    y = NULL,
    size = "Gene hits",
    color = "-log10(FDR)"
  ) +

  theme_classic(base_size = 14)
  
  ggsave("GOBP_neargenelncRNA.svg")
       write.table(GOBP_neargenelncRNA, "GOBP_neargenelncRNA.txt" )
                            
###############################################################################

# gencode/ensembl suele tener columna gene_type o biotype
mRNA_counts <- counts %>%
  filter(gene_type == "protein_coding")  # solo genes codificantes

rownames(mRNA_counts) <- mRNA_counts$Geneid
mRNA_counts <- mRNA_counts[, grep("_sort.bam", colnames(mRNA_counts))]





# ejemplo: correlación expresión lncRNA vs genes
cor_mat <- cor(t(expr_lnc), t(expr_mRNA), method = "spearman")

# seleccionar genes asociados
targets <- names(which(cor_mat["Xist", ] > 0.4))

# enrichment clásico
library(clusterProfiler)
enrich <- enrichGO(
  gene = targets,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP"
)


#####################################################################################################################################
# Script Graphos Known lncRNA


Functional pathway annotation of Known lncRNAs
Differentially expressed lncRNAs with described function were analyzed independently for the Contralateral and Ipsilateral conditions. For each lncRNA, associated biological pathways and predicted gene targets were obtained from functional enrichment analyses. The resulting annotation table included lncRNA identifiers, pathway terms, and gene targets.
To reduce redundancy and improve interpretability, pathway terms were consolidated into higher-order functional categories using case-insensitive keyword matching implemented in R (v4.x). Pathways containing “Wnt” were classified as WNT; those containing “cytoskeleton” or “actin cytoskeleton” as Cytoskeleton-related; “NF-kB” or “Nfkb” as NFkB; “EMT” as EMT; “apoptosis” or “p53” as Apoptosis; “imprinting” as Imprinting; “inflam” as Inflammation; and “neuronal” or “neurogenesis” as Neurogenesis. Pathways not matching predefined categories were assigned to Other.
Functional classification was performed separately for Contralateral and Ipsilateral datasets using custom scripts in R (tidyverse framework).
Network construction
Networks were constructed independently for Contralateral and  Ipsilateral conditions using shared functional classification and lncRNA–target relationships. Network objects were generated using the igraph package (v1.x).
Two types of edges were defined: 1. lncRNA–target edges, generated by reshaping the annotation table to connect each lncRNA to its associated pathway terms and predicted gene targets. 2.Functional similarity edges, generated between pairs of lncRNAs belonging to the same functional category. Pairwise combinations were computed using a self-join strategy within each functional class. Self-interactions were excluded, and duplicate edges were removed by retaining only lexicographically ordered pairs (lncRNA₁ < lncRNA₂), ensuring an undirected network structure. Edges corresponding to the Other category were excluded to minimize non-informative associations.
All edges were merged into a unified edge table containing source node, target node, and edge type. Nodes were defined as the unique set of lncRNAs and targets present in the final edge list. No artificial pathway nodes were introduced; therefore, all nodes correspond to biological entities.
Network visualization and modular organization
Networks were visualized using a force-directed layout (Fruchterman–Reingold algorithm) implemented in ggraph (v2.x). Edge color was mapped to edge type to distinguish functional similarity relationships from lncRNA–target associations. Edge width was adjusted to emphasize functional similarity edges relative to lncRNA–target edges. Nodes were colored according to node type (lncRNA or target gene).
Contralateral and  Ipsilateral networks were constructed and visualized independently to enable condition-specific comparison of functional architecture and module organization.



library(tidyverse)
library(igraph)
library(ggraph)

functional classification
lnc_path_ipsis <- lnc_ipsis %>%
  dplyr::select(lncRNA, pathway) %>%
  dplyr::mutate(functional_class = case_when(
    lncRNA %in% c("H19") | grepl("Wnt", pathway, ignore.case = TRUE) ~ "WNT",
    lncRNA %in% c("MEXIS", "ARHGAP20OS", "CELRR") | grepl(""cytoskeleton|actin cytoskeleton", pathway, ignore.case = TRUE) ~ "Cytoskeleton-related",
    grepl("NF-kB|Nfkb", pathway, ignore.case = TRUE) ~ "NFkB",
    grepl("EMT", pathway, ignore.case = TRUE) ~ "EMT",
    grepl("apoptosis|p53", pathway, ignore.case = TRUE) ~ "Apoptosis",
    grepl("imprinting", pathway, ignore.case = TRUE) ~ "Imprinting",
    grepl("inflam", pathway, ignore.case = TRUE) ~ "Inflammation",
    grepl("neuronal|neurog", pathway, ignore.case = TRUE) ~ "Neurogenesis",
    #grepl("cytoskeleton|actin cytoskeleton", pathway, ignore.case = TRUE) ~ "Cytoskeleton-related",
    TRUE ~ "Other"
  ))

# ---Edges lncRNA -> target ---
edges_ipsis <- lnc_ipsis %>%
  pivot_longer(cols = c(pathway, gene_target),
               names_to = "type",
               values_to = "target") %>%
  dplyr::select(from = lncRNA, to = target) %>%
  dplyr::mutate(type = "lnc_target")

# ---Edges funcionales entre lncRNA ---
pairs_ipsis <- lnc_path_ipsis %>%
  dplyr::inner_join(lnc_path_ipsis, by = "functional_class", relationship = "many-to-many") %>%
  dplyr::filter(lncRNA.x < lncRNA.y) %>%
  dplyr::filter(functional_class != "Other") %>%
  dplyr::select(from = lncRNA.x, to = lncRNA.y, functional_class) %>%
  dplyr::distinct() %>%
  dplyr::mutate(type = functional_class)

# ---Conexión especial Cytoskeleton-related ---
cyto_edges <- tibble(
  from = c("MEXIS", "CELRR"),
  to   = c("Arhgap20", "Arhgap20"),
  type = "Cytoskeleton-related"
)

# ---Combinar todos los edges ---
edges_ipsis_all <- bind_rows(
  edges_ipsis,
  pairs_ipsis %>% dplyr::select(from, to, type),
  cyto_edges
)

# ---Crear nodos ---
nodes_ipsis <- tibble(
  name = unique(c(edges_ipsis_all$from, edges_ipsis_all$to))
) %>%
  dplyr::mutate(node_type = ifelse(name %in% lnc_ipsis$lncRNA, "lncRNA", "target"))

# ---Grafo ---
graph_ipsis <- igraph::graph_from_data_frame(edges_ipsis_all, vertices = nodes_ipsis)

#################################################################################################################

References Databases

                                                                                   link
1                              http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/
2                                              http://jaspar.genereg.net/html/DOWNLOAD/
3                                                                                      
4                                        http://amp.pharm.mssm.edu/lib/cheadownload.jsp
5                                                      http://www.ncbi.nlm.nih.gov/geo/
6                                          http://genome.ucsc.edu/ENCODE/downloads.html
7                                   https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways
8                               http://www.wikipathways.org/index.php/Download_Pathways
9                                                      http://www.ncbi.nlm.nih.gov/geo/
10                                                    http://www.kegg.jp/kegg/download/
11                                                     http://www.ncbi.nlm.nih.gov/geo/
12            http://www.targetscan.org/cgi-bin/targetscan/data_download.cgi?db=vert_61
13                                                        http://amp.pharm.mssm.edu/X2K
14                                            https://pubmed.ncbi.nlm.nih.gov/22110038/
15                             http://software.broadinstitute.org/gsea/msigdb/index.jsp
16                                                         http://biogps.org/downloads/
17                                                         http://biogps.org/downloads/
18                                             http://www.human-phenotype-ontology.org/
19                                                   http://www.roadmapepigenomics.org/
20                                     http://amp.pharm.mssm.edu/lib/keacommandline.jsp
21                                                https://www.nursa.org/nursa/index.jsf
22                                  http://mips.helmholtz-muenchen.de/genre/proj/corum/
23                                     http://amp.pharm.mssm.edu/lib/keacommandline.jsp
24                                                  http://www.broadinstitute.org/cmap/
25                                                  http://www.broadinstitute.org/cmap/
26                                                        http://www.omim.org/downloads
27                                                        http://www.omim.org/downloads
28                                            http://mint.bio.uniroma2.it/download.html
29                            http://www.broadinstitute.org/gsea/msigdb/collections.jsp
30                            http://www.broadinstitute.org/gsea/msigdb/collections.jsp
31                                                     http://www.ncbi.nlm.nih.gov/geo/
32                                                     http://www.ncbi.nlm.nih.gov/geo/
33                                                     http://www.ncbi.nlm.nih.gov/geo/
34                                       https://portals.broadinstitute.org/ccle/home\n
35                                                         http://biogps.org/downloads/
36                                                        https://www.proteomicsdb.org/
37                                            http://www.humanproteomemap.org/index.php
38                                                         http://www.hmdb.ca/downloads
39                                          ftp://ftp.ebi.ac.uk/pub/databases/interpro/
40                                                            http://www.brain-map.org/
41                                         http://genome.ucsc.edu/ENCODE/downloads.html
42                                         http://genome.ucsc.edu/ENCODE/downloads.html
43                                                      http://www.koehn.embl.de/depod/
44                                                            http://www.brain-map.org/
45                                         http://genome.ucsc.edu/ENCODE/downloads.html
46                                               http://www.broadinstitute.org/achilles
47                                               http://www.broadinstitute.org/achilles
48                                  https://cgap.nci.nih.gov/Pathways/BioCarta_Pathways
49                                                                 http://humancyc.org/
50                                                    http://www.kegg.jp/kegg/download/
51                                                            http://www.pantherdb.org/
52                              http://www.wikipathways.org/index.php/Download_Pathways
53                                                     http://www.maayanlab.net/ESCAPE/
54                                               http://www.ncbi.nlm.nih.gov/homologene
55                                                     http://www.ncbi.nlm.nih.gov/geo/
56                                                     http://www.ncbi.nlm.nih.gov/geo/
57                                                     http://www.ncbi.nlm.nih.gov/geo/
58                                              https://grants.nih.gov/grants/oer.htm\n
59                                                     http://www.ncbi.nlm.nih.gov/geo/
60                                                    http://amp.pharm.mssm.edu/Enrichr
61                                                     http://www.ncbi.nlm.nih.gov/geo/
62                                                     http://www.ncbi.nlm.nih.gov/geo/
63                                                    http://amp.pharm.mssm.edu/Enrichr
64                                                      http://www.ncbi.nlm.nih.gov/gap
65                                                                     https://clue.io/
66                                                                     https://clue.io/
67                                                           http://www.gtexportal.org/
68                                                           http://www.gtexportal.org/
69                                                     http://www.ncbi.nlm.nih.gov/geo/
70                                                     http://www.ncbi.nlm.nih.gov/geo/
71                                                     http://www.ncbi.nlm.nih.gov/geo/
72                                                     http://www.ncbi.nlm.nih.gov/geo/
73                                                     http://www.ncbi.nlm.nih.gov/geo/
74                                                     http://www.ncbi.nlm.nih.gov/geo/
75                                                     http://www.ncbi.nlm.nih.gov/geo/
76                                                     http://www.ncbi.nlm.nih.gov/geo/
77                                                                     https://clue.io/
78                                                                     https://clue.io/
79                                                                     https://clue.io/
80                                                                     https://clue.io/
81                                                    http://www.kegg.jp/kegg/download/
82                              http://www.wikipathways.org/index.php/Download_Pathways
83                                                                                     
84                                                     http://www.ncbi.nlm.nih.gov/geo/
85                                                     http://www.ncbi.nlm.nih.gov/geo/
86                                   http://cgap.nci.nih.gov/Pathways/BioCarta_Pathways
87                                                                 http://humancyc.org/
88  https://www.ndexbio.org/index.html#/networkset/7bc65b82-2a2f-11ed-ac45-0ac135e8bacf
89                                                    http://www.pantherdb.org/pathway/
90                                                https://ntp.niehs.nih.gov/drugmatrix/
91                                                    http://amp.pharm.mssm.edu/Enrichr
92                                                         http://proteincomplexes.org/
93                                                        http://tissues.jensenlab.org/
94                                                     http://www.ncbi.nlm.nih.gov/geo/
95                                                   http://compartments.jensenlab.org/
96                                                       http://diseases.jensenlab.org/
97                                                      http://bioplex.hms.harvard.edu/
98                                                     http://amp.pharm.mssm.edu/archs4
99                                                     http://amp.pharm.mssm.edu/archs4
100                                                    http://amp.pharm.mssm.edu/archs4
101                                                    http://amp.pharm.mssm.edu/archs4
102                                                    http://amp.pharm.mssm.edu/archs4
103                                                         http://sys-myo.rhcloud.com/
104                                                  http://mirtarbase.mbc.nctu.edu.tw/
105                                                          http://www.targetscan.org/
106                                                   http://amp.pharm.mssm.edu/Enrichr
107                                                   http://amp.pharm.mssm.edu/Enrichr
108                                                   http://amp.pharm.mssm.edu/Enrichr
109                                       http://tanlab.ucdenver.edu/DSigDB/DSigDBv1.0/
110                                                    http://www.ncbi.nlm.nih.gov/geo/
111                                       http://hgdownload.cse.ucsc.edu/downloads.html
112                                                https://amp.pharm.mssm.edu/geneshot/
113                                     https://www.ncbi.nlm.nih.gov/gene/about-generif
114                                     https://www.ncbi.nlm.nih.gov/gene/about-generif
115                                                https://amp.pharm.mssm.edu/geneshot/
116                                                      http://www.subcellbarcode.org/
117                                                          https://www.ebi.ac.uk/gwas
118                                                       https://www.wikipathways.org/
119                                                       https://www.wikipathways.org/
120                                                    https://www.grnpedia.org/trrust/
121                                                                https://www.kegg.jp/
122                                                                https://www.kegg.jp/
123                                                     https://www.ebi.ac.uk/interpro/
124                                                              https://pfam.xfam.org/
125                                                                 https://depmap.org/
126                                                                 https://depmap.org/
127                                               https://www.ukbiobank.ac.uk/tag/gwas/
128                                                   https://tripod.nih.gov/bioplanet/
129                                               https://www.ncbi.nlm.nih.gov/clinvar/
130                                                        http://pheweb.sph.umich.edu/
131                                                            https://www.disgenet.org
132                                            http://lincs.hms.harvard.edu/kinomescan/
133                                             https://portals.broadinstitute.org/ccle
134                                                       https://www.proteomicsdb.org/
135                                                  https://amp.pharm.mssm.edu/lnchub/
136                                                                http://phipster.org/
137                                           http://www.transgene.ru/disease-pathways/
138                                                                                    
139                                                  https://amp.pharm.mssm.edu/covid19
140                            https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp 
141                                                     https://maayanlab.cloud/Enrichr
142                                               https://toxico.nibiohn.go.jp/english/
143                                                       https://portal.brain-map.org/
144     https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/
145                                                                https://www.kegg.jp/
146                                                       https://www.wikipathways.org/
147                               https://hubmapconsortium.github.io/ccf-asct-reporter/
148                                                        http://www.geneontology.org/
149                                                        http://www.geneontology.org/
150                                                        http://www.geneontology.org/
151                                                     http://www.informatics.jax.org/
152                                               http://biocc.hrbmu.edu.cn/CellMarker/
153                                                           http://www.orphadata.org/
154                                                    https://maayanlab.cloud/covid19/
155                                                               https://panglaodb.se/
156                                               https://azimuth.hubmapconsortium.org/
157                                            https://www.ncbi.nlm.nih.gov/gap/phegeni
158                                                             https://gtexportal.org/
159                                                             https://www.hdinhd.org/
160                                                             https://www.hdinhd.org/
161                               https://hubmapconsortium.github.io/ccf-asct-reporter/
162                                                      https://fantom.gsc.riken.jp/6/
163                          https://github.com/nybell/drugsets/tree/main/DATA/GENESETS
164                                      https://tabula-sapiens-portal.ds.czbiohub.org/
165                                                      https://maayanlab.cloud/chea3/
166                   https://appyters.maayanlab.cloud/#/Gene_Expression_T2D_Signatures
167                                     https://maayanlab.cloud/sigcom-lincs/#/Download
168                                     https://maayanlab.cloud/sigcom-lincs/#/Download
169                                               https://tabula-muris.ds.czbiohub.org/
170                                                  https://reactome.org/download-data
171                                                        https://www.syngoportal.org/
172                                                             https://www.glygen.org/
173                                                            https://drugcentral.org/
174                                                     https://www.mousephenotype.org/
175                                              https://www.metabolomicsworkbench.org/
176                                  https://www.nature.com/articles/s41587-022-01539-0
177                                         https://kinase-library.phosphosite.org/site
178                                                        https://gtexportal.org/home/
179                                                        http://www.geneontology.org/
180                                                        http://www.geneontology.org/
181                                                        http://www.geneontology.org/
182                                                     https://pfocr.wikipathways.org/
183                                                          https://www.ebi.ac.uk/gwas
184                                                    http://gedipnet.bicnirrh.res.in/
185                                             https://magnet-winterlab.herokuapp.com/
186                                               https://azimuth.hubmapconsortium.org/
187                                                              https://rummagene.com/
188                                                              https://rummagene.com/
189                                                              https://rummagene.com/
190                                                           https://motrpac-data.org/
191                                                       https://www.wikipathways.org/
192                                                                 https://depmap.org/
193                                                        https://www.syngoportal.org/
194                                         http://bio-bigdata.hrbmu.edu.cn/CellMarker/
195                                                                  https://dgidb.org/
196                                                     http://www.informatics.jax.org/
197                                                       https://www.wikipathways.org/
198                                                       https://www.wikipathways.org/
199                                         https://kinase-library.phosphosite.org/site
200                                                https://perturbatlas.kratoss.site/#/
201                                                  https://reactome.org/download-data
202                                                        http://www.geneontology.org/
203                                                        http://www.geneontology.org/
204                                                        http://www.geneontology.org/
205                                                                                    
206                                                                                    
207                                                      http://diseases.jensenlab.org/
208                                                      http://diseases.jensenlab.org/
209                                           https://compartments.jensenlab.org/Search
210                                                https://perturbatlas.kratoss.site/#/
211                                                            https://gwps.wi.mit.edu/
212                                                https://tissues.jensenlab.org/Search
213                                           https://compartments.jensenlab.org/Search
214                                                            https://gwps.wi.mit.edu/
215                                                https://tissues.jensenlab.org/Search
216                                                               https://rummageo.com/
217                                                               https://rummageo.com/
218                                                          https://www.ebi.ac.uk/gwas
219                                               https://www.ncbi.nlm.nih.gov/clinvar/
220                                                                  https://cm4ai.org/
221                                                           https://jaspar.elixir.no/
222                                                           https://jaspar.elixir.no/
223                                                           https://jaspar.elixir.no/
224                                                          https://carcinogenome.org/
225                                                                https://www.kegg.jp/
