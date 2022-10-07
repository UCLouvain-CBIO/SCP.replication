
R_HOME  = "$(shell R RHOME)"
R = "$(R_HOME)/bin/R"

website_shell:
	make init
	make home
	make reference
	make articles_index
	make news

init:
	${R} --quiet -e "pkgdown::init_site()";

home:
	${R} --quiet -e "pkgdown::build_home()";

reference:
	${R} --quiet -e "pkgdown::build_reference()";

articles_index:
	${R} --quiet -e "pkgdown::build_articles_index()";

news:
	${R} --quiet -e "pkgdown::build_news()";
