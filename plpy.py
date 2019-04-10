class plpy:
# substitute for postgres' plpy when using outside postgres

    def notice(self, msg):
        print msg;
