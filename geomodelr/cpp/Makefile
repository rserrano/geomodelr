
cpp: cpp.so
	cp cpp.so ..

cpp.so: basic.o section.o match.o limiter.o model.o geomodel.o faults.o polygon.o
	g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wl,-Bsymbolic-functions -Wl,-z,relro -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security  basic.o section.o match.o limiter.o model.o geomodel.o faults.o polygon.o -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/ -lboost_python -o cpp.so

polygon.o: polygon.cpp polygon.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c polygon.cpp -std=c++14
basic.o: basic.cpp basic.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c basic.cpp -std=c++14
section.o: section.cpp section.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c section.cpp -std=c++14
match.o: match.cpp match.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c match.cpp -std=c++14
limiter.o: limiter.cpp limiter.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c limiter.cpp -std=c++14
model.o: model.cpp model.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c model.cpp -std=c++14
geomodel.o: geomodel.cpp geomodel.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c geomodel.cpp -std=c++14
faults.o: faults.cpp faults.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c faults.cpp -std=c++14

clean:
	rm *.o
	rm cpp.so
	rm ../cpp.so
