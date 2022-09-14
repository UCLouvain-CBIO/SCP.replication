
R_HOME  = "$(shell R RHOME)"
R = "$(R_HOME)/bin/R"

docker_build:
	docker build -t cvanderaa/scp_replication_docker:v1 .

website_update:
	make site
	make docs_push

site:
	make init
	make home
	make reference
	make articles
	make news

init:
	${R} -e "pkgdown::init_site()";

home:
	${R} -e "pkgdown::build_home()";

reference:
	${R} -e "pkgdown::build_reference()";

articles:
    ## Compile vignettes using cvanderaa/scp_replication_docker:v1
	# docker pull cvanderaa/scp_replication_docker:v1
	docker run -e PASSWORD=scp \
	           -v /home/cvanderaa/PhD/SCP.replication/docs:/home/rstudio/SCP.replication/docs \
			   -it cvanderaa/scp_replication_docker:v1 \
			   make articles_in_v1
	${R} -e "pkgdown::build_articles_index()";
	
articles_in_v1:
	${R} -e "pkgdown::build_article('SCoPE2')";
	# ${R} -e "pkgdown::build_article('zhu2019EL')";
	${R} -e "pkgdown::build_article('schoof2021')";
	${R} -e "pkgdown::build_article('liang2020')";
	${R} -e "pkgdown::build_article('williams2020_tmt')";
	# ${R} -e "pkgdown::build_article('leduc2022')";
	# ${R} -e "pkgdown::build_article('derks2022')";

news:
	${R} -e "pkgdown::build_news()";

docs_push:
	git add docs/*
	git com7mit -m "Compiled vignettes"
	git push