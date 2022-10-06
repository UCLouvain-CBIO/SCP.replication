
R_HOME  = "$(shell R RHOME)"
R = "$(R_HOME)/bin/R"

website:
	make init
	make home
	make reference
	make articles
	make news

init:
	${R} --quiet -e "pkgdown::init_site()";

home:
	${R} --quiet -e "pkgdown::build_home()";

reference:
	${R} --quiet -e "pkgdown::build_reference()";

articles:
    ## Compile vignettes using cvanderaa/scp_replication_docker:v1
	docker pull cvanderaa/scp_replication_docker:v1
	make articles_in_v1
	## Uncomment below to insert new version
	## docker pull cvanderaa/scp_replication_docker:v2
	## make articles_in_v2 
	${R} -e "pkgdown::build_articles_index()";
	
articles_in_v1:
	docker run -e PASSWORD=scp \
		-v `pwd`:/home/rstudio/SCP.replication/ \
		-it cvanderaa/scp_replication_docker:v1 \
		R --quiet -e "setwd('SCP.replication');pkgdown::build_article('leduc2022');pkgdown::build_article('brunner2022');pkgdown::build_article('derks2022');pkgdown::build_article('liang2020');pkgdown::build_article('schoof2021');pkgdown::build_article('SCoPE2');pkgdown::build_article('williams2020_tmt');pkgdown::build_article('zhu2019EL')"

## Uncomment to insert new version
## articles_in_v2:
##  docker run -e PASSWORD=scp \
##		-v `pwd`:/home/rstudio/SCP.replication/ \
##		-it cvanderaa/scp_replication_docker:v2 \
##      R --quiet -e "setwd('SCP.replication');pkgdown::build_article('NEW_VIGNETTE')\""

news:
	${R} --quiet -e "pkgdown::build_news()";
