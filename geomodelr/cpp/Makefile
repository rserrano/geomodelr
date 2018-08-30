
cpp: cpp.so
	cp cpp.so ..

cpp.so: basic.o section.o match.o model.o geomodel.o faults.o speed_up.o polygon.o isosurfaces_vdb.o feflow.o
	g++ -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -Wl,-Bsymbolic-functions -Wl,-z,relro -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security  basic.o section.o match.o model.o geomodel.o faults.o speed_up.o polygon.o isosurfaces_vdb.o feflow.o -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu/ -lboost_python -lopenvdb -lgmp -lmpfr -lCGAL -o cpp.so

polygon.o: polygon.cpp polygon.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c polygon.cpp -std=c++11
basic.o: basic.cpp basic.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c basic.cpp -std=c++11
section.o: section.cpp section.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c section.cpp -std=c++11
match.o: match.cpp match.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c match.cpp -std=c++11
model.o: model.cpp model.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c model.cpp -std=c++11
geomodel.o: geomodel.cpp geomodel.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c geomodel.cpp -std=c++11
faults.o: faults.cpp faults.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c faults.cpp -std=c++11
speed_up.o: speed_up.cpp speed_up.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c speed_up.cpp -std=c++11
isosurfaces_vdb.o: isosurfaces_vdb.cpp isosurfaces_vdb.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c isosurfaces_vdb.cpp -std=c++11
feflow.o: feflow.cpp feflow.hpp
	g++ -pthread -DNDEBUG -g -fwrapv -O2 -Wall -fno-strict-aliasing -Wdate-time -D_FORTIFY_SOURCE=2 -g -fstack-protector-strong -Wformat -Werror=format-security -fPIC -I/usr/local/include -I/usr/include/python2.7 -c feflow.cpp -std=c++11