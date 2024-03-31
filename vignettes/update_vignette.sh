cp ~/research_code/sceptre-book/ondisc.qmd ~/research_code/ondisc/vignettes/ondisc.Rmd

sed -i -e '/::: callout-note/,/:::/d' ondisc.Rmd
sed -i -e 's|\[frequently asked questions page\](faq.qmd)|\[frequently asked questions page\](https://timothy-barry.github.io/sceptre-book/faq.html)|g' ondisc.Rmd
sed -i -e 's/{#[^}]*}//g' ondisc.Rmd
sed -i -e 's/[^.]*@sec-[^.]*\.//g' ondisc.Rmd
echo "
\`\`\`{r}
library(sessioninfo); session_info()
\`\`\`" >> ondisc.Rmd
sed -i '' '1i\
--- \
title: "Getting started with ondisc" \
output: rmarkdown::html_vignette \
vignette: > \
  %\\VignetteIndexEntry{Getting started with ondisc} \
  %\\VignetteEngine{knitr::rmarkdown} \
  %\\VignetteEncoding{UTF-8} \
--- \
' ondisc.Rmd

rm ondisc.Rmd-e
