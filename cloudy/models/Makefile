%.out: %.in
	/Users/will/Work/CLOUDY/gitlab-cloudy/source/cloudy.exe -p $* 

OUTFILES := $(patsubst %.in,%.out,$(wildcard w3-*.in))
TEST_OUTFILES := $(patsubst %.in,%.out,$(wildcard test*.in))

all: $(OUTFILES)

test: $(TEST_OUTFILES)
