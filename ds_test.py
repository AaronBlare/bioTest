import os
import sys
import logging
from openai import OpenAI
from openai import APIError, APIConnectionError, RateLimitError, AuthenticationError

LOG_FILE = "errors.log"
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.ERROR,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

API_KEY = ''
if not API_KEY:
    print("Ошибка: не задана переменная окружения DEEPSEEK_API_KEY", file=sys.stderr)
    sys.exit(1)

client = OpenAI(
    api_key=API_KEY,
    base_url="https://api.deepseek.com/" 
)

MODEL = "deepseek-chat"

if len(sys.argv) > 1:
    user_message = " ".join(sys.argv[1:])
else:
    user_message = input("Введите ваш вопрос: ").strip()
    if not user_message:
        print("Пустой вопрос. Завершение работы.", file=sys.stderr)
        sys.exit(0)

try:
    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "user", "content": user_message}
        ],
        stream=False,
        timeout=30.0
    )
    
    answer = response.choices[0].message.content
    print("\nОтвет DeepSeek:\n", answer)

except AuthenticationError as e:
    error_msg = f"Ошибка аутентификации: неверный API ключ. {str(e)}"
    logging.error(error_msg)
    print(f"Ошибка: {error_msg}", file=sys.stderr)
    sys.exit(1)

except RateLimitError as e:
    error_msg = f"Превышен лимит запросов. {str(e)}"
    logging.error(error_msg)
    print(f"Ошибка: {error_msg}", file=sys.stderr)
    sys.exit(1)

except APIConnectionError as e:
    error_msg = f"Проблема соединения с API DeepSeek. {str(e)}"
    logging.error(error_msg)
    print(f"Ошибка: {error_msg}", file=sys.stderr)
    sys.exit(1)

except APIError as e:
    error_msg = f"Ошибка API DeepSeek (статус, внутренняя ошибка и т.п.): {str(e)}"
    logging.error(error_msg)
    print(f"Ошибка: {error_msg}", file=sys.stderr)
    sys.exit(1)

except Exception as e:
    error_msg = f"Неизвестная ошибка: {str(e)}"
    logging.error(error_msg)
    print(f"Ошибка: {error_msg}", file=sys.stderr)
    sys.exit(1)