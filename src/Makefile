ARCH=$(shell uname)
ifeq ($(ARCH),Darwin)
	include Makefile.osx
else ifeq ($(ARCH),Linux)
	include Makefile.linux
else
	$(warning Unsupported architecture. Trying Makefile.linux)
	include Makefile.linux
endif
