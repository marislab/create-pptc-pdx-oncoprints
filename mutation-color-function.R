#library(gplots)
#library(ComplexHeatmap)

col2hex(c('mediumseagreen', 'orange', 'green2', 'hotpink',
          'lightslateblue', 'midnightblue', 'gray33', 'turquoise3',
          'lightpink', 'lightsalmon', 'khaki1',
          'purple', 'black', 'dodgerblue2', 'red', 'cornflowerblue', 'lightcoral', 'red', 'mistyrose', 'lightsteelblue1'))

col = c("Missense_Mutation" = "#009E73", 
        "Nonsense_Mutation" = "#000000",
        "Frame_Shift_Del" = "#56B4E9", 
        "Frame_Shift_Ins" = "#CC79A7", 
        "Splice_Site" = "#F0E442",
        "Translation_Start_Site" = "#191970",
        "Nonstop_Mutation" = "#545454",
        "In_Frame_Del" = "#CAE1FF",
        "In_Frame_Ins" = "#FFE4E1",
        "Stop_Codon_Ins" = "#CC79A7",
        "Start_Codon_Del" = "#56B4E9",
        "Fusion" = "#A020F0",
        "Multi_Hit" = "#E69F00",
        #"Deletion" = "#0072B2",
        "Hom_Deletion" = "#0072B2",
        "Hem_Deletion" = "#9FB6CD",
        "Amplification" = "#D55E00",
        "Loss" = "#0072B2",
        "Gain" = "#D55E00",
        "High_Level_Gain" = "#FF0000",
        "Multi_Hit_Fusion" = "381373")

mut.labels = c("Missense Mutation", 
        "Nonsense Mutation",
        "Frameshift Deletion", 
        "Frameshift Insertion", 
        "Splice Site Mutation",
        "Translation Start Site",
        "Nonstop Mutation",
        "Inframe Deletion",
        "Inframe Insertion",
        "Stop Codon Insertion",
        "Start Codon Deletion",
        "RNA Fusion",
        "Multi-Hit",
        "Homozygous Deletion",
        "Hemizygous Deletion",
        "Focal Amplification",
        "Loss",
        "Gain",
        "High Level Gain",
        "Multi-Hit RNA Fusion")


alter_fun = function(x, y, w, h, v) {
            if(v["Missense_Mutation"]) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["Missense_Mutation"], col = bg.col))
            if(v["Nonsense_Mutation"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Nonsense_Mutation"], col = bg.col))
            if(v["Frame_Shift_Del"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Frame_Shift_Del"], col = bg.col))
            if(v["Frame_Shift_Ins"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Frame_Shift_Ins"], col = bg.col))
            if(v["Splice_Site"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Splice_Site"], col = bg.col))
            if(v["Translation_Start_Site"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Translation_Start_Site"], col = bg.col))
            if(v["Nonstop_Mutation"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Nonstop_Mutation"], col = bg.col))
            if(v["In_Frame_Del"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["In_Frame_Del"], col = bg.col))
            if(v["In_Frame_Ins"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["In_Frame_Ins"], col = bg.col))
            if(v["Stop_Codon_Ins"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Stop_Codon_Ins"], col = bg.col))
            if(v["Start_Codon_Del"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Start_Codon_Del"], col = bg.col))
            if(v["Fusion"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Fusion"], col = bg.col))
            if(v["Multi_Hit"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Multi_Hit"], col = bg.col))
            if(v["Hom_Deletion"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Hom_Deletion"], col = bg.col))
            if(v["Hem_Deletion"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Hem_Deletion"], col = bg.col))
            if(v["Amplification"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Amplification"], col = bg.col))
            if(v["Loss"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Loss"], col = bg.col))
            if(v["Gain"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Gain"], col = bg.col))
            if(v["High_Level_Gain"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["High_Level_Gain"], col = bg.col))
             if(v["Multi_Hit_Fusion"]) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["Multi_Hit_Fusion"], col = bg.col))
}
