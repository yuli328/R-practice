#来源：《R数据科学》
#Time:2026.03.02
#PKGs:tidyverse,含ggplot2\tibble\readr\purrr\dplyr
install.packages("tidyverse")
library(tidyverse)
#基础概念：数据框mpg、图形属性aes()、映射标度图层、分面（每一组都用单独一张小图展示）、几何对象（散点图、折线图、箱线图等）、统计变换、
#位置调整geom_bar的position参数设置为“identity”、“fill”、“dodge”，示例：
# 创建数据：3个班级，2种性别，人数
df <- data.frame(
  class = c("一班", "一班", "二班", "二班", "三班", "三班"),
  gender = c("男", "女", "男", "女", "男", "女"),
  count = c(20, 18, 25, 22, 15, 17)
)
#position = "stack"
ggplot(df, aes(x = class, y = count, fill = gender)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "position = 'dodge' | 并列对比")
#position = "fill" 
ggplot(df, aes(x = class, y = count, fill = gender)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "position = 'fill' | 百分比占比") +
  scale_y_continuous(labels = scales::percent) # 显示百分比
#position = "identity"
ggplot(df, aes(x = class, y = count, fill = gender)) +
  geom_bar(stat = "identity", position = "identity", alpha = 0.5) +
  labs(title = "position = 'identity' | 重叠（加透明度才能看清）")

#总结：绘图关键点：数据框解析，厘清数据框变量，确定合适的几何对象；映射为何种图形属性，是否要分面；标度和统计变化要明确。

#练习：

#一、
a<- mpg
ggplot(data=a)
?mpg  #查询变量
ggplot(data=mpg)+
  geom_point(mapping=aes(x=cyl,y=hwy)) #气缸数和油耗

#二、区别
ggplot(data=mpg)+
  geom_point(mapping=aes(x=displ,y=hwy,color="blue")) 
ggplot(data=mpg)+
  geom_point(mapping=aes(x=displ,y=hwy),color="blue") 
str(a)
#连续变量/分类变量映射为 color、size 和 shape，混乱/易区分
?geom_point

#三、分面
facet_grid(drv ~ cyl)
ggplot(data=mpg)+
  geom_point(mapping=aes(x=drv,y=cyl))
?facet_grid #facet_grid()函数分面，参数为公式，~前代表行，后代表列。如果将行来分面就是var~.;如果作为列来分面就是.~var。如果只有一个变量用来分面，也可以用facet_wrap()。
?facet_wrap
ggplot(data=mpg)+
  geom_point(mapping=aes(x=displ,y=hwy))+
  facet_grid(drv~.)
ggplot(data=mpg)+
  geom_point(mapping=aes(x=displ,y=hwy))+
  facet_grid(.~cyl)

#四、几何对象
#折线图：geom_line()箱线图：geom_boxplot()直方图：geom_histogram()分区图：facet_grid()\facet_wrap()
ggplot(
  data=mpg,
  mapping=aes(x=displ,y=hwy,color=drv)
  )+geom_point()+
  geom_smooth(se=FALSE)#全局映射 #se 是否显示误差范围
#点图，反映x与y的关系，根据drv的不同显示不同颜色的点。有三条平滑的回归线。

ggplot(data=mpg,mapping=aes(x=displ,y=hwy))+geom_point()+
  geom_smooth()

ggplot()+
  geom_point(
    data=mpg,
    mapping=aes(x=displ,y=hwy)
  )+
  geom_smooth(
    data=mpg,
    mapping=aes(x=displ,y=hwy)
  )
#为什么展示的图是一样的？**仅是映射的图层不同

#思考编图
ggplot(data=mpg,mapping=aes(x=displ,y=hwy))+
  geom_point()+
  geom_smooth()

ggplot(data=mpg,mapping=aes(x=displ,y=hwy,group=drv))+
  geom_point()+
  geom_smooth(se=FALSE)

ggplot(data=mpg,mapping=aes(x=displ,y=hwy,color=drv))+
  geom_point()+
  geom_smooth(se=FALSE)

ggplot(data=mpg,mapping=aes(x=displ,y=hwy))+
  geom_point(mapping=aes(color=drv))+
  geom_smooth(se=FALSE)

ggplot(data=mpg,mapping=aes(x=displ,y=hwy))+
  geom_point(mapping=aes(color=drv))+
  geom_smooth(aes(linetype=drv),se=FALSE)

ggplot(data=mpg,
       mapping=aes(x=displ,y=hwy,color=drv))+
  geom_point()

#五、统计变换
?stat_summary()
ggplot(data=diamonds)+
  geom_pointrange(
    aes(x=cut,y=depth),
    stat="summary",
    fun.min=min,
    fun.max=max,
    fun=median
  )
?geom_col() #条形图来显示数据的值
?geom_bar() #条形图来显示数据的比例

#比例条形图为什么需要设定 group = 1？
ggplot(data=diamonds)+
  geom_bar(mapping=aes(x=cut,y=..prop..))
ggplot(data=diamonds)+
  geom_bar(mapping=aes(x=cut,fill=color,y=prop..)
  )
#分类变量所有类别为一个组，然后分别计算每个类别的比例

#六、位置调整
?geom_jitter() #直接显示可能被覆盖的点
?geom_count()  #黑点大小来显示此处事件密度
ggplot(mpg,aes(x=drv,y=displ))+geom_boxplot()




