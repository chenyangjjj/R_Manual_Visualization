
graphics.off()
rm(list=ls())

# 加载所需要的包
library(ggpubr)


#加载数据集ToothGrowth
data("ToothGrowth")
df1 <- ToothGrowth
head(df1)


p = ggboxplot(df1, x="dose", y="len", color = "dose", 
               palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
               add = "jitter",shape = "dose")#增加了jitter点，点shape由dose映射
p1 = p + stat_compare_means(method = "anova")



## 指定组别的比较
my_comparisons <- list( c("0.5", "1"), c("0.5", "2"), c("1", "2") )
p = ggboxplot(df1, x = "dose", y = "len",
                     color = "dose",add = "jitter", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
p1 = p + stat_compare_means(method = "t.test",comparisons = my_comparisons)

ggsave("C:/Users/18543/Desktop/R_组会/p2.tiff",p2, units="in", width=6.5, height=5)


## 指定参考组的比较
p = ggboxplot(df1, x = "dose", y = "len",
          color = "dose",add = "jitter", palette = c("#00AFBB", "#E7B800", "#FC4E07"))
p3 = p +  stat_compare_means(method = "anova") + 
  stat_compare_means(label = "p.signif", method = "t.test",ref.group = "0.5")   

ggsave("C:/Users/18543/Desktop/R_组会/p3.tiff",p3, units="in", width=6.5, height=5)



Plot_Comparison_Baseline_ROI_FTP_SUVR_Stages_CSF_PET_Abeta = 
  ggarrange(p3,p3,p3,p3,
    nrow=2,ncol = 2,
    font.label = list(size = 14, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend="bottom")





## 定义不同的比较方法
# stat_compare_means(method = "wilcox.test","t.test","anova",……)


## 其他图形的绘制
x1 = ggviolin(df1, x = "dose", y = "len", fill = "dose",
                      palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                      add = "boxplot", add.params = list(fill = "white")) + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
               stat_compare_means(label.y = 50)    

ggbarplot(df, x = "dose", y = "len", fill = "dose",
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             add = "boxplot", add.params = list(fill = "white"))
ggsave("C:/Users/18543/Desktop/R_组会/x1.tiff",x1, units="in", width=6.5, height=5)


##Deviation graphs
data("mtcars")
dfm <- mtcars
# Convert the cyl variable to a factor
dfm$cyl <- as.factor(dfm$cyl)
# Add the name colums
dfm$name <- rownames(dfm)
# Inspect the data
head(dfm[, c("name", "wt", "mpg", "cyl")])


# Inspect the data
ggbarplot(dfm, x = "name", y = "mpg",
          fill = "cyl",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in dscending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90)          # Rotate vertically x axis texts



dfm$mpg_z <- (dfm$mpg -mean(dfm$mpg))/sd(dfm$mpg)
dfm$mpg_grp <- factor(ifelse(dfm$mpg_z < 0, "low", "high"), 
                      levels = c("low", "high"))

ggbarplot(dfm, x = "name", y = "mpg_z",
          fill = "mpg_grp",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "MPG z-score",
          legend.title = "MPG Group",
          rotate = TRUE,
          ggtheme = theme_minimal())




