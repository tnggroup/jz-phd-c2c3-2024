#Standardised Mean Difference template
#Source: https://handbook-5-1.cochrane.org/chapter_9/9_2_3_2_the_standardized_mean_difference.htm


#Let's say you have a dataframe with your individual level data - toy example
df.cases <- data.frame(
  case=1,
  detachment=rnorm(n = 200, mean= 3.5, sd = 1),
  negative_affect=rnorm(n = 200, mean= 4.5, sd = 1),
  antagonism=rnorm(n = 200, mean= 2.5, sd = 1),
  disinhibition=rnorm(n = 200, mean= 3.5, sd = 1),
  psychoticism=rnorm(n = 200, mean= 3, sd = 1.5)
)

df.controls <- data.frame(
  case=0,
  detachment=rnorm(n = 200, mean= 3, sd = 1.5),
  negative_affect=rnorm(n = 200, mean= 4, sd = 1),
  antagonism=rnorm(n = 200, mean= 2.5, sd = 1),
  disinhibition=rnorm(n = 200, mean= 3, sd = 1.5),
  psychoticism=rnorm(n = 200, mean= 2, sd = 1.5)
  )

df <- rbind(df.cases,df.controls) #df now contains both 200 cases and 200 controls


#standardised mean difference in detachment between cases and controls
## standardised mean difference is the mean difference in measurement relative to the overall spread (standard deviation) of that measurement - different measurements have different spread

detachment.smd <- mean(df[df$case==1,]$detachment-df[df$case==0,]$detachment)/sd(df$detachment)
##Uses the pooled standard deviation of detachment across both groups!
detachment.smd

negative_affect.smd <- mean(df[df$case==1,]$negative_affect-df[df$case==0,]$negative_affect)/sd(df$negative_affect)
##Uses the pooled standard deviation of negative affect across both groups!
negative_affect.smd

#Now the comparison between SMD in detachment and SMD in negative affect, both relative measurements to their respective spread, is more informative

