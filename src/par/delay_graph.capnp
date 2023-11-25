@0xb1182a5115fb9537;

struct Fanin {
    regCnt @0 :UInt32;
    combDelay @1 :UInt32;
    driverId @2 :UInt64;
}

struct Vertex{
    id @0 :UInt64;
    componentId @1 :UInt64;
    outportId @2 :UInt16;
    fanins @3 :List(Fanin);
}

struct Graph{
    vertices @0 :List(Vertex);
    primaryInputId @1 :List(UInt64);
    primaryOutputId @2 :List(UInt64);
}