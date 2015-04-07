
import smtplib

from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

mailto_list=["nanyang.xu@gmail.com"]

mail_host="202.38.64.8"
default_subject = 'computation result notify'
default_title = 'your task is done, and this is the brieve result'

def send_result(to_list,result,mail_user, mail_pass,imagefile=None):
    me=mail_user
    msg = MIMEMultipart()
    msg['Subject'] = default_subject
    msg['From'] = me
    msg['To'] = ";".join(to_list)
    msg.attach(MIMEText(result))
    if imagefile != None:
        fp = open(imagefile,'rb')
        msg.attach(MIMEImage(fp.read()))
        fp.close()
        
    try:
        s = smtplib.SMTP()
        s.connect(mail_host)
        s.login(mail_user,mail_pass)
        s.sendmail(me, to_list, msg.as_string())
        s.close()
        return True
    except Exception, e:
        print str(e)
        return False
    
if __name__ == '__main__':
    if send_result(mailto_list,"subject","content"):
        print "Sent OK"
    else:
        print "Sent failure"

