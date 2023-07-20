function status = myoutfcn(t,app,flag,tstar,fname)
  switch flag
    case 'init'
      status = 0;
    case ''
      if (t(end) == tstar)
        app = app(:,end);
        save([fname,'.mat'],'app')
      end
      status = 0;
    case 'done'
      status = 1;
  end
end
