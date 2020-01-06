import sphere_database : Database;

void main()
{
    auto database = new Database( "databases/B407_M2.nc");
    database.eval( 1000., 45.2, 75.0, 0.15, -6.0, 0.004);
    database.eval( 1000., -180.0, 180.0, 0.15, -6.0, 0.004);
}
