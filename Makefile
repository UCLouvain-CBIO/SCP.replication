
R_HOME  = "$(shell R RHOME)"
R = "$(R_HOME)/bin/R"

docker_build:
	docker build -t cvanderaa/scp_replication_docker:v1 .

website_update:
	make website
	make docs_push

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
	# docker pull cvanderaa/scp_replication_docker:v1
	make articles_in_v1
	${R} -e "pkgdown::build_articles_index()";
	
articles_in_v1:
	docker run -e PASSWORD=scp \
	           -v /home/cvanderaa/PhD/SCP.replication/docs:/home/rstudio/SCP.replication/docs \
			   -it cvanderaa/scp_replication_docker:v1 \
			   R --quiet -e "pkgdown::build_article('SCoPE2'); pkgdown::build_article('zhu2019EL'); pkgdown::build_article('schoof2021'); pkgdown::build_article('liang2020'); pkgdown::build_article('williams2020_tmt')"
	# ${R} -e "pkgdown::build_article('leduc2022')";
	# ${R} -e "pkgdown::build_article('derks2022')";

news:
	${R} --quiet -e "pkgdown::build_news()";

docs_push:
	git add docs/*
	git commit -m "Compiled vignettes"
	git push