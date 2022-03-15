from datetime import datetime


def make_monthly_chunks(start: datetime, end: datetime) -> "list[dict]":
    chunks = []
    
    for year in range(start.year, end.year):
        
        for month in range(1, 13):
            
            chunk = {'year': str(year), 
                     'month': "{:02d}".format(month),
                     'time': ["{:02d}:00".format(H) for H in range(0, 24)]}
            if year == start.year and month < start.month:
                continue
            
            elif year == start.year and month == start.month:
                day = ["{:02d}".format(d) for d in range(start.day, 32)]
            
            elif year == end.year and month == end.month:
                day = ["{:02d}".format(d) for d in range(1, end.day + 1)]
            
            elif year == end.year and month > end.month:
                continue
            
            else: 
                day = ["{:02d}".format(d) for d in range(1, 32)]
            
            chunk['day'] = day

            if len(chunk['day'])  == 1:
                chunk['day'] = chunk['day'][0]
            
            chunks.append(chunk)
            
    return chunks